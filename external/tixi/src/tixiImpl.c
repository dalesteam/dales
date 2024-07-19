/*
* Copyright (C) 2015 German Aerospace Center (DLR/SC)
*
* Created: 2010-08-13 Markus Litz <Markus.Litz@dlr.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
/**
 * @file   tixiImpl.c
 * @author Hans-Peter Kersken <Hans-Peter.Kersken@dlr.de>
 * @date   Mon Mar 14 16:19:56 2005
 * @brief  Implementation of the TIXI-C-API
 * @todo add not an Element to tixiGetTextAttribute
 */

/*
* Altered by Sven Werchner (KIT) in order to create a minimal version for ICON-ART
* Note: When XML-functionality for ICON-ART is expanding, 
*       this TIXI-Version probably needs to be (re-)expanded as well
*/


#include <assert.h>

#include "tixiInternal.h"
#include "tixiUtils.h"
#include "tixi.h"
#include "tixi_version.h"

//static xmlNsPtr nameSpace = NULL;

void tixiDefaultMessageHandler(MessageType type, const char* msg);

/**
  Pointer to the head of the list of documents managed by TIXI
*/
TixiDocumentListEntry *documentListHead = NULL;

TixiPrintMsgFnc tixiMessageHandler = tixiDefaultMessageHandler;

/**
  Route all messages to the internal message handler
 */
void printMsg(MessageType type, const char* message, ...)
{
  char buffer[2048];
  char* extra = NULL;
  int len = 0;

  va_list varArgs;
  va_start(varArgs, message);
  if ((len = vsnprintf(buffer, 2048, message, varArgs)) >= 2048) {
    // message is longer than 2048 bytes,
    // we must allocate
    extra = (char*) malloc((len+1) * sizeof(char));

    vsnprintf(extra, len+1, message, varArgs);
    if (tixiMessageHandler) {
      tixiMessageHandler(type, extra);
    }

    if (extra) {
      free(extra);
    }
  }
  else {
    if (tixiMessageHandler) {
      tixiMessageHandler(type, buffer);
    }
  }

  va_end(varArgs);
}

int _initialized = 0;

void xmlErrorHandler(void * ctx, const char *message, ...) {
  char buffer[2048];
  char* extra = NULL;
  int len = 0;

  va_list varArgs;
  va_start(varArgs, message);
  if ((len = vsnprintf(buffer, 2048, message, varArgs)) >= 2048) {
    // message is longer than 2048 bytes,
    // we must allocate
    extra = (char*) malloc((len+1) * sizeof(char));

    vsnprintf(extra, len+1, message, varArgs);
    tixiMessageHandler(MESSAGETYPE_ERROR, extra);

    if (extra) {
      free(extra);
    }
  }
  else {
    tixiMessageHandler(MESSAGETYPE_ERROR, buffer);
  }

  va_end(varArgs);
}

static void tixiInit(void)
{
  if (!_initialized) {
    printMsg(MESSAGETYPE_STATUS, "TiXI initialized\n");
    xmlSetGenericErrorFunc(NULL, xmlErrorHandler);
    _initialized = 1;
  }
}


/**
  gives the tixi version number
*/
DLL_EXPORT char* tixiGetVersion()
{
  static char version[] = TIXI_VERSION_STRING;
  return version;
}



/**
 *  Opens the file and sets up the TixiDocument datastructure.
 */
DLL_EXPORT ReturnCode tixiOpenDocumentRecursive(const char *xmlFilename, TixiDocumentHandle *handle, OpenMode oMode)
{
  /* this opens the XML-file and checks if it is well formed */

  TixiDocument *document = NULL;
  xmlDocPtr xmlDocument = NULL;
  FILE *file = NULL;
  ReturnCode returnValue = -1;

  tixiInit();
  checkLibxml2Version();

  xmlKeepBlanksDefault(0);
  xmlIndentTreeOutput = 1;

  assert(xmlFilename);

  file = fopen(xmlFilename, "r");
  if (!file) {
    printMsg(MESSAGETYPE_ERROR, "Error: Unable to open file \"%s\".\n", xmlFilename);
    return OPEN_FAILED;
  }
  else {
    fclose(file);
  }

  xmlDocument = xmlReadFile(xmlFilename, NULL, 0);

  if (xmlDocument) {

    document = (TixiDocument *) malloc(sizeof(TixiDocument));
    document->xmlFilename = (char *) malloc(strlen(xmlFilename) * sizeof(char) + 1);
    strcpy(document->xmlFilename, xmlFilename);
    strip_dirname(xmlFilename, &document->dirname, &document->filename);
    document->validationFilename = NULL;
    document->docPtr = xmlDocument;
    document->currentNode = NULL;
    document->isValid = UNDEFINED;
    document->status = OPENED;
    document->memoryListHead = NULL;
    document->memoryListTail = NULL;
    document->hasIncludedExternalFiles = 1;
    document->usePrettyPrint = 1;
    document->uidListHead = NULL;
    addDocumentToList(document, &(document->handle));
    *handle = document->handle;
    returnValue = SUCCESS; /*?*/

    if (oMode == OPENMODE_RECURSIVE) {
      int count = 0;
      document->hasIncludedExternalFiles = 1;

      returnValue = openExternalFiles(document, &count);
      if (returnValue != SUCCESS){
        printMsg(MESSAGETYPE_ERROR, "Error %d in including external files into tixiDoument.\n", returnValue);
        removeDocumentFromList(*handle);
        freeTixiDocument(document);
        document = NULL;
      }
    }
  }
  else {
    printMsg(MESSAGETYPE_ERROR, "Error: \"%s\" is not a wellformed XML-file.\n", xmlFilename);
    returnValue = NOT_WELL_FORMED;
  }
  return returnValue;
}


DLL_EXPORT ReturnCode tixiOpenDocument(const char *xmlFilename, TixiDocumentHandle *handle)
{
  return tixiOpenDocumentRecursive(xmlFilename, handle, OPENMODE_PLAIN);
}



DLL_EXPORT ReturnCode tixiCloseDocument(TixiDocumentHandle handle)
{
  TixiDocument *document = getDocument(handle);

  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle in tixiCloseDocument.\n");
    return INVALID_HANDLE;
  }

  if (removeDocumentFromList(handle) == FAILED) {
    if (document) {
      freeTixiDocument(document);
      document = NULL;
    }
    return CLOSE_FAILED;
  }

  if (document) {
    freeTixiDocument(document);
    document = NULL;
  }

  return SUCCESS;
}



DLL_EXPORT ReturnCode tixiGetTextElement(const TixiDocumentHandle handle, const char *elementPath, char **text)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  xmlNodePtr element = NULL;
  ReturnCode error = SUCCESS;

  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }

  xmlDocument = document->docPtr;

  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);

  if (!error) {

    xmlNodePtr children = element->children;
    char *textPtr = NULL;

    textPtr = (char *) xmlNodeListGetString(xmlDocument, children, 0);
    if ( textPtr ) {
      *text = (char *) malloc((strlen(textPtr) + 1) * sizeof(char));
      strcpy(*text, textPtr);
      xmlFree(textPtr);
    } else {
      *text = (char *) malloc(sizeof(char));
      strcpy(*text, "");
    }
    error = addToMemoryList(document, (void *) *text);
    if (!error) {
      xmlXPathFreeObject(xpathObject);
    }
  }

  return error;
}

DLL_EXPORT ReturnCode tixiGetIntegerElement(const TixiDocumentHandle handle, const char *elementPath, int *number)
{
  char *text;
  ReturnCode error = 0;

  error = tixiGetTextElement(handle, elementPath, &text);

  if (error) {
    printMsg(MESSAGETYPE_STATUS, "Error: tixiGetTextElement returns %d in tixiGetIntegerElement.\n", error);
    return error;
  }

  *number = atoi(text);

  return SUCCESS;
}

DLL_EXPORT ReturnCode tixiGetDoubleElement(const TixiDocumentHandle handle, const char *elementPath, double *number)
{
  char *text = NULL;
  ReturnCode error = 0;

  error = tixiGetTextElement(handle, elementPath, &text);

  if (error) {
    printMsg(MESSAGETYPE_STATUS, "Error: tixiGetTextElement returns %d in tixiGetDoubleElement.\n", error);
    return error;
  }

  *number = atof(text);
  return SUCCESS;

  /* use   strtod(nptr, (char **)NULL); to check for errors */
}


DLL_EXPORT ReturnCode tixiGetTextAttribute(const TixiDocumentHandle handle, const char *elementPath,
                                           const char *attributeName, char **text)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  char *textPtr;
  xmlNodePtr element = NULL;
  ReturnCode error = SUCCESS;


  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }

  xmlDocument = document->docPtr;

  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);
  if (!error) {

    textPtr = (char *) xmlGetProp(element, (xmlChar *) attributeName);

    if (textPtr) {
      *text = (char *) malloc((strlen(textPtr) + 1) * sizeof(char));
      strcpy(*text, textPtr);
      xmlFree(textPtr);
      error = addToMemoryList(document, (void *) *text);
      xmlXPathFreeObject(xpathObject);
      return error;
    }
    else {
      /*printMsg(MESSAGETYPE_ERROR,
              "Error: Attribute \"%s\" of element \"%s\" not found.  \n",
              attributeName, elementPath); */
      xmlXPathFreeObject(xpathObject);
      return ATTRIBUTE_NOT_FOUND;
    }
  }

  return error;
}

DLL_EXPORT ReturnCode tixiGetDoubleAttribute(const TixiDocumentHandle handle,
                                             const char *elementPath, const char *attributeName, double *number)
{
  char *text;
  ReturnCode error = 0;



  error = tixiGetTextAttribute(handle, elementPath, attributeName, &text);

  if (error) {
    printMsg(MESSAGETYPE_ERROR, "Error: tixiGetTextAttribute returns %d in tixiGetDoubleAttribute.\n", error);
    return error;
  }

  *number = atof(text);
  /* use   strtod(nptr, (char **)NULL); to check for errors */
  return SUCCESS;
}

DLL_EXPORT ReturnCode tixiGetIntegerAttribute(const TixiDocumentHandle handle,
                                              const char *elementPath, const char *attributeName, int *number)
{
  char *text;
  ReturnCode error = 0;



  error = tixiGetTextAttribute(handle, elementPath, attributeName, &text);

  if (error) {
    printMsg(MESSAGETYPE_STATUS, "Error: tixiGetTextAttribute returns %d in tixiGetIntegerAttribute.\n", error);
    return error;
  }

  *number = atoi(text);
  /* use   strtod(nptr, (char **)NULL); to check for errors */
  return SUCCESS;
}


DLL_EXPORT ReturnCode tixiGetNamedChildrenCount(const TixiDocumentHandle handle,
                                                const char *elementPath, const char *childName, int *count)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathContextPtr xpathContext = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  xmlNodeSetPtr nodes = NULL;
  int iNode;
  char *childElementPath =
      (char *) malloc(sizeof(char) * (strlen(elementPath) + strlen(childName) + 2));
  char *allChildren = (char *) malloc(sizeof(char) * (strlen(elementPath) + 3));



  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    free(childElementPath);
    free(allChildren);
    return INVALID_HANDLE;
  }

  childElementPath[0] = '\0';
  allChildren[0] = '\0';

  strcat(childElementPath, elementPath);
  strcat(childElementPath, "/");
  strcat(childElementPath, childName);

  strcat(allChildren, elementPath);
  strcat(allChildren, "/*");

  *count = 0;

  xmlDocument = document->docPtr;

  xpathContext = xmlXPathNewContext(xmlDocument);

  if (!xpathContext) {
    printMsg(MESSAGETYPE_ERROR, "Error: unable to create new XPath context\n");
    free(childElementPath);
    free(allChildren);
    return FAILED;
  }

  /* first check parent */
  xpathObject = xmlXPathEvalExpression((xmlChar *) elementPath, xpathContext);

  if (!xpathObject) {
    printMsg(MESSAGETYPE_ERROR, "Error: unable to evaluate xpath expression \"%s\"\n", elementPath);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return INVALID_XPATH;
  }

  if (xmlXPathNodeSetIsEmpty(xpathObject->nodesetval)) {
    xmlXPathFreeObject(xpathObject);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return ELEMENT_NOT_FOUND;
  }

  nodes = xpathObject->nodesetval;
  assert(nodes);

  if (nodes->nodeNr > 1) {
    printMsg(MESSAGETYPE_ERROR,
             "Error: Element chosen by XPath \"%s\" expression is not unique. \n", elementPath);
    free(childElementPath);
    free(allChildren);
    xmlXPathFreeObject(xpathObject);
    xmlXPathFreeContext(xpathContext);
    return ELEMENT_PATH_NOT_UNIQUE;
  }

  xmlXPathFreeObject(xpathObject);

  /* check if there are children at all */


  xpathObject = xmlXPathEvalExpression((xmlChar *) allChildren, xpathContext);

  if (!xpathObject) {
    printMsg(MESSAGETYPE_ERROR, "Error: unable to evaluate xpath expression \"%s\"\n", allChildren);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return INVALID_XPATH;
  }


  if (xmlXPathNodeSetIsEmpty(xpathObject->nodesetval)) {

    /* parent has no child at all, return child count 0 */
    *count = 0;
    xmlXPathFreeObject(xpathObject);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return SUCCESS;
  }

  xmlXPathFreeObject(xpathObject);

  /* now check child */

  xpathObject = xmlXPathEvalExpression((xmlChar *) childElementPath, xpathContext);

  if (!xpathObject) {
    printMsg(MESSAGETYPE_ERROR, "Error: unable to evaluate xpath expression \"%s\"\n", childElementPath);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return INVALID_XPATH;
  }

  if (xmlXPathNodeSetIsEmpty(xpathObject->nodesetval)) {

    /* parent has no child with name childName , return child count 0 */
    *count = 0;
    xmlXPathFreeObject(xpathObject);
    xmlXPathFreeContext(xpathContext);
    free(childElementPath);
    free(allChildren);
    return SUCCESS;
  }

  nodes = xpathObject->nodesetval;
  assert(nodes);


  for (iNode = 0; iNode < nodes->nodeNr; iNode++) {

    assert(nodes->nodeTab[iNode]);

    if (nodes->nodeTab[iNode]->type == XML_ELEMENT_NODE) {

      if (!strcmp(childName, (char *) nodes->nodeTab[iNode]->name)) {
        (*count)++;
      }
    }
  }

  xmlXPathFreeContext(xpathContext);
  xmlXPathFreeObject(xpathObject);
  free(childElementPath);
  free(allChildren);
  return SUCCESS;
}




DLL_EXPORT ReturnCode   tixiGetChildNodeName(const TixiDocumentHandle handle, const char *elementPath,  int index, char **text)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlNodePtr element = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  int error = SUCCESS;

  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }

  xmlDocument = document->docPtr;

  if(index <= 0){
    return INDEX_OUT_OF_RANGE;
  }

  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);
  xmlXPathFreeObject(xpathObject);

  if(!error){
    xmlNodePtr child = element->children;
    int pos = 1;

    while(child && pos < index){
      // Ignore DTD nodes
      if (child->type != XML_DTD_NODE) {
        pos++;
      }
      child = child->next;
    }

    if(pos != index || !child){
      return INDEX_OUT_OF_RANGE;
    }

    // return node value according to dom specification: http://www.w3schools.com/dom/dom_nodetype.asp
    if(child->type == XML_TEXT_NODE){
      *text = (char *) malloc(10 * sizeof(char));
      strcpy(*text, "#text");
    }
    else if(child->type == XML_CDATA_SECTION_NODE){
      *text = (char *) malloc(20 * sizeof(char));
      strcpy(*text, "#cdata-section");
    }
    else if(child->type == XML_COMMENT_NODE){
      *text = (char *) malloc(10 * sizeof(char));
      strcpy(*text, "#comment");
    }
    else {
      // get name
      *text = (char *) malloc((strlen((char*)child->name) + 3) * sizeof(char));
      strcpy(*text,  (char*) child->name);
    }
    error = addToMemoryList(document, (void *) *text);
  }
  return error;
}

DLL_EXPORT ReturnCode tixiGetNumberOfChilds(const TixiDocumentHandle handle, const char *elementPath, int* nChilds)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  xmlNodePtr element = NULL;
  ReturnCode error = SUCCESS;


  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }
  xmlDocument = document->docPtr;


  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);
  xmlXPathFreeObject(xpathObject);

  if (!error) {
    xmlNodePtr children = element->children;
    *nChilds = 0;
    while(children){
      // Ignore DTD nodes, we don't select them with xpath
      if (children->type != XML_DTD_NODE) {
        (*nChilds)++;
      }
      children = children->next;
    }
  }

  return error;
}


void tixiDefaultMessageHandler(MessageType type, const char *message)
{
  // only show errors and warnings by default
  if (type < MESSAGETYPE_STATUS) {
    fputs(message, stderr);
  }
}

DLL_EXPORT ReturnCode tixiGetNumberOfAttributes(const TixiDocumentHandle handle, const char *elementPath, int* nAttributes)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  xmlNodePtr element = NULL;
  ReturnCode error = SUCCESS;


  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }
  xmlDocument = document->docPtr;

  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);
  xmlXPathFreeObject(xpathObject);

  if (!error) {
    xmlAttrPtr attr = element->properties;
    *nAttributes = 0;

    while(attr){
      attr = attr->next;
      (*nAttributes)++;
    }
  }

  return error;
}

DLL_EXPORT ReturnCode tixiGetAttributeName(const TixiDocumentHandle handle, const char *elementPath, int attrIndex, char** attrName)
{
  TixiDocument *document = getDocument(handle);
  xmlDocPtr xmlDocument = NULL;
  xmlXPathObjectPtr xpathObject = NULL;
  xmlNodePtr element = NULL;
  ReturnCode error = SUCCESS;


  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle.\n");
    return INVALID_HANDLE;
  }
  xmlDocument = document->docPtr;

  if(attrIndex <= 0){
    return INDEX_OUT_OF_RANGE;
  }

  error = checkElement(xmlDocument, elementPath, &element, &xpathObject);
  xmlXPathFreeObject(xpathObject);


  if (!error) {
    xmlAttrPtr attr = element->properties;
    int pos = 1;

    while(attr && pos < attrIndex){
      attr = attr->next;
      pos++;
    }

    if(pos != attrIndex || !attr){
      return INDEX_OUT_OF_RANGE;
    }

    // get name
    *attrName = (char *) malloc((strlen((char*)attr->name) + 3) * sizeof(char));
    strcpy(*attrName,  (char*)attr->name);
    error = addToMemoryList(document, (void *) *attrName);
    }
   
  return error;
}
   

DLL_EXPORT ReturnCode tixiGetDocumentPath(TixiDocumentHandle handle, char** documentPath)
{
  TixiDocument *document = NULL;

  if (!documentPath) {
    printMsg(MESSAGETYPE_ERROR, "Error: Null Pointer in tixiGetDocumentPath.\n");
    return FAILED;
  }

  document = getDocument(handle);
  if (!document) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid document handle in tixiGetDocumentPath.\n");
    return INVALID_HANDLE;
  }

  *documentPath = document->xmlFilename;

  return SUCCESS;
}

