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

#ifndef TIXI_INTERNAL_H
#define TIXI_INTERNAL_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file   tixiInternal.h
 * @author Hans-Peter Kersken <Hans-Peter.Kersken@dlr.de>
 * @date   Tue Mar 15 12:06:33 2005
 * @brief  Header file for routines used to implement the API.
 * This file cotains declaration for auxillary routines used to implement the tixi API.
 */
#include "libxml/xpath.h"
#include "tixi.h"
#include "tixiData.h"



typedef enum {
  SUCESS,
  MEMORY_ALLOCATION_FAILED
} InternalReturnCode;


/**
 * @brief Adds documents to list of managed documents.
 *
 *
 * @param document pointer to the document to be added to the list of managed documents
 *
 * @return handle to the document
 */
InternalReturnCode addDocumentToList(TixiDocument* document, TixiDocumentHandle* handle);

/**
 * @brief Removes document from the list of managed documents.
 *
 *
 * @param document pointer to the document to be added to the list of managed documents
 *
 * @return error code
 */
ReturnCode removeDocumentFromList(TixiDocumentHandle handle);

/**
 * @brief Retrives a document.
 *
 *
 * @param (in) handle
 *
 * @return handle to the document or NULL if handle not in the list
 */
DLL_EXPORT TixiDocument* getDocument(TixiDocumentHandle handle);

/**
 @brief Adds pointer to memory allcoated by libxml to list of managed memory.

 @param (in) memory pointer to the memory location
 */
InternalReturnCode addToMemoryList(TixiDocument* document, void* memory);

/**
  @brief Frees all memory locations pointed to by the pointers in the list
         and removes all list entries.

  @param document (in) a pointer to a TixiDocument structure
 */
InternalReturnCode clearMemoryList(TixiDocument* document);

/**
  @brief Frees the memory used by the document

  @param document (in) a pointer to a TixiDocument structure
*/
void freeTixiDocument(TixiDocument* document);


/**
  @brief Checks if the given element path is valid

  @param xmlDocument (in) pointer to an libxml2 document
  @param elementPath (in) path to the element to be check
  @param element (out) pointer to the XML-node pointed to be element path
  @param xpathObject (out) libxml2 internal object pointer. Has to be freed by the calling routine.
  @return
    - SUCCESS element exists and is unique
    - FAILED internal error
    - INVALID_XPATH
    - ELEMENT_NOT_FOUND
    - NOT_AN_ELEMENT
 */
DLL_EXPORT ReturnCode checkElement(const xmlDocPtr xmlDocument, const char* elementPath, xmlNodePtr* element,
                                   xmlXPathObjectPtr* xpathObject);

/**
  brief Check if libxml2 version used to build TIXI.
*/
void checkLibxml2Version();

/**
  @brief builds a string using snprintf to avoid buffer overflows. The resulting buffer
         has to be freed by the caller.
*/
DLL_EXPORT char* buildString(const char* format, ...);

/**
  @brief Open external xml files and merge them into the tree.

  External xml files could be integrated into the xml tree by linking them into the main xml file.
  The user now has only one big xml file.
  For example:

    <wings>
        <airfoils>
           <externaldata>
                <path>./WingSections/</path>
                <filename>VFW614-W-1.xml</filename>
                <filename>VFW614-W-2.xml</filename>
                ...
            </externaldata>
            <airfoil>
                <name>VFW614 Seitenleitwerksprofil</name>
                <coordinates>
                    <point><x>1.0000000</x><y>0.0000000</y><z>0.0000000</z></point>
                      <point><x>0.9795687</x><y>0.0023701</y><z>0.0000000</z></point>
                      ...
                </coordinates>
            </airfoil>
        </airfoils>
    </wings>

  @param TixiDocument (in) a TIXI document with a opened main-xml file.
  @param int (out) The number of files that where includet.

  @return
    - SUCCESS
    - FAILED internal error
    - OPEN_FAILED
*/
ReturnCode openExternalFiles(TixiDocument* aTixiDocument, int* number);


#ifdef __cplusplus
}
#endif

#endif /* TIXI_INTERNAL_H */
