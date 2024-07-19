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

/*
* Altered by Sven Werchner (KIT) in order to create a minimal version for ICON-ART
* Note: When XML-functionality for ICON-ART is expanding, 
*       this TIXI-Version probably needs to be (re-)expanded as well
*/

#include <assert.h>

#include "tixi.h"
#include "tixiData.h"
#include "tixiInternal.h"
#include "xpathFunctions.h"
#include "uidHelper.h"
#include "tixiUtils.h"

/**
   @file Auxiliary routines used to implement the interface.
*/



extern void printMsg(MessageType type, const char* message, ...);

InternalReturnCode clearMemoryList(TixiDocument* document)
{

  TixiMemoryListEntry* current = document->memoryListHead;

  while (current) {

    TixiMemoryListEntry* next = (TixiMemoryListEntry*) current->next;
    free(current->memory);
    free(current);
    current = next;
  }
  document->memoryListHead = NULL;
  document->memoryListTail = NULL;
  return SUCCESS;
}

void freeTixiDocument(TixiDocument* document)
{
  if (document->xmlFilename) {
    free(document->xmlFilename);
    document->xmlFilename = NULL;
  }

  if (document->validationFilename) {
    free(document->validationFilename);
    document->validationFilename = NULL;
  }

  if (document->dirname) {
    free(document->dirname);
    document->dirname = NULL;
  }

  if (document->filename) {
    free(document->filename);
    document->filename = NULL;
  }
  clearMemoryList(document);
  uid_clearUIDList(document);

  xmlFreeDoc(document->docPtr);

  free(document);
}

InternalReturnCode addDocumentToList(TixiDocument* document, TixiDocumentHandle* handle)
{

  TixiDocumentListEntry* currentEntry = documentListHead;
  static int handleCounter = 0;

  if (documentListHead) {

    while (currentEntry->next) {
      currentEntry = (TixiDocumentListEntry*) currentEntry->next;
    }

    currentEntry->next = (struct TixiDocumentListEntry*)
        malloc(sizeof(TixiDocumentListEntry));

    if (!currentEntry->next) {

      return MEMORY_ALLOCATION_FAILED;
    }

    currentEntry = (TixiDocumentListEntry*) currentEntry->next;
  }
  else {

    currentEntry = (TixiDocumentListEntry*) malloc(sizeof(TixiDocumentListEntry));

    if (currentEntry) {
      documentListHead = currentEntry;
    }
    else {
      return MEMORY_ALLOCATION_FAILED;
    }
  }

  currentEntry->document = document;
  currentEntry->next = NULL;

  handleCounter++;
  *handle = handleCounter;

  return SUCCESS;
}


ReturnCode removeDocumentFromList(TixiDocumentHandle handle)
{

  TixiDocumentListEntry* currentEntry = documentListHead;
  TixiDocumentListEntry* previousEntry = NULL;

  ReturnCode returnValue = SUCCESS;

  if (documentListHead) {

    while (currentEntry) {

      if (currentEntry->document->handle != handle) {
        previousEntry = currentEntry;
        currentEntry = (TixiDocumentListEntry*) currentEntry->next;
      }
      else {
        if (previousEntry) {
          previousEntry->next = currentEntry->next;
          free(currentEntry);
          currentEntry = NULL;
          returnValue = SUCCESS;
        }
        else {
          documentListHead = (TixiDocumentListEntry*) currentEntry->next;
          free(currentEntry);
          currentEntry = NULL;
          returnValue = SUCCESS;
        }
      }
    }
  }
  else {
    returnValue = FAILED;
  }

  return returnValue;
}

TixiDocument* getDocument(TixiDocumentHandle handle)
{

  TixiDocumentListEntry* currentEntry = documentListHead;
  TixiDocument* returnValue = NULL;

  if (currentEntry) {
    while (currentEntry && returnValue == NULL) {

      if (currentEntry->document->handle != handle) {
        currentEntry = (TixiDocumentListEntry*) currentEntry->next;
      }
      else {
        returnValue = currentEntry->document;
      }
    }
  }
  else {
    returnValue = NULL;
  }

  return returnValue;
}

InternalReturnCode addToMemoryList(TixiDocument* document, void* memory)
{

  TixiMemoryListEntry* currentEntry = document->memoryListTail;

  if (currentEntry) {

    currentEntry->next = (struct TixiMemoryListEntry*) malloc(sizeof(TixiMemoryListEntry));

    if (currentEntry->next) {
      currentEntry = (TixiMemoryListEntry*) currentEntry->next;
      document->memoryListTail = currentEntry;
    }
    else {
      return MEMORY_ALLOCATION_FAILED;
    }
  }
  else {

    currentEntry = (TixiMemoryListEntry*) malloc(sizeof(TixiMemoryListEntry));

    if (currentEntry) {
      document->memoryListHead = currentEntry;
      document->memoryListTail = currentEntry;
    }
    else {
      return MEMORY_ALLOCATION_FAILED;
    }
  }

  currentEntry->memory = memory;
  currentEntry->next = NULL;

  return SUCCESS;
}


ReturnCode checkElement(const xmlDocPtr xmlDocument, const char* elementPathDirty,
                        xmlNodePtr* element, xmlXPathObjectPtr* xpathObject)
{

  xmlXPathContextPtr xpathContext = NULL;
  xmlNodeSetPtr nodes = NULL;
  char elementPath[1024];

  *xpathObject = NULL;

  /* remove trailing slash */
  strncpy(elementPath, elementPathDirty, 1024);
  if(strlen(elementPath) > 1 && elementPath[strlen(elementPath)-1] == '/') {
    elementPath[strlen(elementPath)-1] = '\0';
  }

  /* Create xpath evaluation context */
  xpathContext = xmlXPathNewContext(xmlDocument);
  if (!xpathContext) {
    printMsg(MESSAGETYPE_ERROR, "Error: Unable to create new XPath context.\n");
    return FAILED;
  }

  /* Evaluate Expression */
  *xpathObject = xmlXPathEvalExpression((xmlChar*) elementPath, xpathContext);
  if (!(*xpathObject)) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid XPath expression \"%s\"\n", elementPath);
    xmlXPathFreeContext(xpathContext);
    return INVALID_XPATH;
  }

  if (xmlXPathNodeSetIsEmpty((*xpathObject)->nodesetval)) {
    char * errorStr = buildString("Error: element %s not found!", elementPath);
    xmlXPathFreeContext(xpathContext);
    xmlXPathFreeObject(*xpathObject);
    *xpathObject = NULL;

    printMsg(MESSAGETYPE_STATUS, errorStr);
    free(errorStr);
    return ELEMENT_NOT_FOUND;
  }

  nodes = (*xpathObject)->nodesetval;
  assert(nodes);

  if (nodes->nodeNr > 1) {
    printMsg(MESSAGETYPE_ERROR,
             "Error: Element chosen by XPath \"%s\" expression is not unique. \n", elementPath);
    xmlXPathFreeContext(xpathContext);
    xmlXPathFreeObject(*xpathObject);
    *xpathObject = NULL;
    return ELEMENT_PATH_NOT_UNIQUE;
  }

  assert(nodes->nodeTab[0]);

  if (nodes->nodeTab[0]->type == XML_ELEMENT_NODE || nodes->nodeTab[0]->type == XML_DOCUMENT_NODE) {
    *element = nodes->nodeTab[0];
    xmlXPathFreeContext(xpathContext);
    return SUCCESS;
  }
  else {
    printMsg(MESSAGETYPE_ERROR,
             "Error: XPath expression \"%s\"does not point to an element node.\n", elementPath);
    xmlXPathFreeContext(xpathContext);
    xmlXPathFreeObject(*xpathObject);
    *xpathObject = NULL;
    return NOT_AN_ELEMENT;
  }
}


void checkLibxml2Version()
{

  static int testPerformed = 0;

  if (!testPerformed) {
    LIBXML_TEST_VERSION;
  }

}

char* buildString(const char* format, ...)
{

  int bufferLength = 10;        /* initial guess for the buffer length */
  int nChars = -1;
  va_list variableList;
  char* buffer = NULL;

  buffer = (char*) malloc(bufferLength * sizeof(char));

  if (buffer) {

    while (nChars < 0) {

      /* try to print into buffer */

      va_start(variableList, format);
      nChars = VSNPRINTF(buffer, bufferLength, format, variableList);
      va_end(variableList);

      if (nChars < bufferLength && nChars > -1) {
        break;
      }

      /* msvc and glibc < 2.1: VSNPRINTF returns a negative integer on failure */
      if (nChars < 0) {

        bufferLength *= 2;
        //PRINT_DEBUG2("nChar = %d buffer length = %d, reallocating\n", nChars, bufferLength);
      }
      /*
               glibc 2.0: VSNPRINTF returns the length ( numbers of chars) of the result
               without the trailing \0 if the length would not exceed bufferLength.
            */
      else {

        bufferLength = nChars + 1;
        //PRINT_DEBUG2("nChar = %d buffer length = %d, reallocating\n", nChars, bufferLength);
      }

      nChars = -1;
      buffer = (char*) realloc(buffer, bufferLength * sizeof(char));

      if (!buffer) {
        break;
      }
    }
  }
  return buffer;
}

char* loadExternalFileToString(const char* filename)
{
  if (isURIPath(filename) != 0) {
    // local file
    return loadFileToString(filename);
  }
  else if (string_startsWith(filename, "file://") == 0) {
    char* localPath = uriToLocalPath(filename);
    char* result = NULL;
    if (!localPath) {
      return NULL;
    }
    result = loadFileToString(localPath);

    free(localPath);
    return result;
  }
  else {
    return NULL; //curlGetURLInMemory(filename);
  }
}

ReturnCode openExternalFiles(TixiDocument* aTixiDocument, int* number)
{
  int iNode = 0;
  int handle = aTixiDocument->handle;
  xmlNodePtr cur = NULL;
  ReturnCode error = SUCCESS;

  assert(aTixiDocument != NULL);
  *number = 0;

  while(1) {
    // loop until there are no externaldata nodes included

    xmlXPathObjectPtr xpathObject = XPathEvaluateExpression(aTixiDocument->docPtr, "//externaldata");
    xmlNodeSetPtr nodeset = NULL;
    char* externalDataNodeXPath, *externalDataDirectoryXPath, *externalDataDirectory, *resolvedDirectory;
    int externalFileCount = 0;

    if (!xpathObject) {
      // no more external data, stop
      break;
    }

    nodeset = xpathObject->nodesetval;
    if (!nodeset || nodeset->nodeNr < 1) {
      break;
    }

    // goto the first node that is an element
    for (iNode = 0; iNode < nodeset->nodeNr; ++iNode) {
      cur = nodeset->nodeTab[iNode];
      if (cur->type == XML_ELEMENT_NODE) {
        break; // for loop
      }
    }
    if (iNode == nodeset->nodeNr) {
      // no element node found
      xmlXPathFreeObject(xpathObject);
      break; // while loop
    }

    // found external data node
    xmlXPathFreeObject(xpathObject);

    /* get nodes XPath */
    externalDataNodeXPath = (char*) xmlGetNodePath(cur);


    /* now get the subdirectory */
    externalDataDirectoryXPath = buildString("%s/%s", externalDataNodeXPath, EXTERNAL_DATA_NODE_NAME_PATH);

    error = tixiGetTextElement(handle, externalDataDirectoryXPath, &externalDataDirectory);
    free(externalDataDirectoryXPath);
    if (error) {
      printMsg(MESSAGETYPE_ERROR, "Error: openExternalFiles returns %d. No path defined in externaldata node!\n", error);
      xmlFree(externalDataNodeXPath);
      return OPEN_FAILED;
    }

    // resolv data directory (in case of relative paths)
    resolvedDirectory = resolveDirectory(aTixiDocument->dirname, externalDataDirectory);

    /* now get number and names of all external files */
    tixiGetNamedChildrenCount(handle, externalDataNodeXPath, EXTERNAL_DATA_NODE_NAME_FILENAME, &externalFileCount);
    if (externalFileCount == 0) {
      printMsg(MESSAGETYPE_ERROR, "Error: no filename nodes defined in externalData node.\n");
      xmlFree(externalDataNodeXPath);
      free(resolvedDirectory);
      return OPEN_FAILED;
    }

    for (iNode = 1; iNode <= externalFileCount; iNode++) {
      char* externalFileName, *externalFullFileName, *newDocumentString, *fileNameXPath;
      xmlDocPtr xmlDocument = NULL;

      fileNameXPath = buildString("%s/filename[%d]", externalDataNodeXPath, iNode);

      tixiGetTextElement(handle, fileNameXPath, &externalFileName);
      free(fileNameXPath);

      /* Build complete filename */
      externalFullFileName = buildString("%s%s", resolvedDirectory, externalFileName);

      /* open files */
      newDocumentString = loadExternalFileToString(externalFullFileName);
      if (newDocumentString == NULL) {
        printMsg(MESSAGETYPE_ERROR, "\nError in fetching external file \"%s\".\n", externalFullFileName);
        free(externalFullFileName);
        xmlFree(externalDataNodeXPath);
        free(resolvedDirectory);
        return OPEN_FAILED;
      }

      /* now parse the file to DOM */
      xmlDocument = xmlReadMemory(newDocumentString, (int) strlen(newDocumentString), "urlResource", NULL, 0);
      free(newDocumentString);

      if (xmlDocument) {
        xmlNodePtr rootToInsert = xmlDocGetRootElement(xmlDocument);

        xmlNodePtr parent = cur->parent;
        if (parent) {
          xmlChar* nodePathNew = NULL;
          char* dataURI = localPathToURI(externalDataDirectory);
          xmlNodePtr nodeToInsert = xmlDocCopyNode(rootToInsert, aTixiDocument->docPtr, 1);

          /* add metadata to node, to allow saving external node data */
          xmlSetProp(nodeToInsert, (xmlChar*) EXTERNAL_DATA_XML_ATTR_FILENAME, (xmlChar*) externalFileName);

          /* save the sub-directory */
          xmlSetProp(nodeToInsert, (xmlChar*) EXTERNAL_DATA_XML_ATTR_DIRECTORY, (xmlChar*) dataURI);
          free(dataURI);

          /* save the external data node position */
          nodePathNew = xmlGetNodePath(parent);
          xmlSetProp(nodeToInsert, (xmlChar*) EXTERNAL_DATA_XML_ATTR_NODEPATH, nodePathNew);
          xmlFree(nodePathNew);

          /* replace externalData node with xml file's content */
          xmlReplaceNode(cur, nodeToInsert);

          /* file could be loaded and parsed, increase the counter */
          (*number)++;
        }

        xmlFreeDoc(xmlDocument);
      }
      else {
        printMsg(MESSAGETYPE_WARNING,
                 "Document %s will be ignored. No valid XML document!\n",
                 externalFullFileName);

        /* remove external data node */
        xmlUnlinkNode(cur);
      }
      free(externalFullFileName);
    } /* end for files */

    free(resolvedDirectory);
    free(externalDataNodeXPath);
    xmlFreeNode(cur);
  }


  if (*number == 0) {
    printMsg(MESSAGETYPE_WARNING, "WARNING: Unable to load any externaldata files.\n");
  }

  return SUCCESS;
}
