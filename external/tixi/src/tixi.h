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

#ifdef __cplusplus
extern "C" {
#endif

#if defined(WIN32)
#if defined (TIXI_EXPORTS)
#define DLL_EXPORT __declspec (dllexport)
#else
#define DLL_EXPORT
#endif
#else
#define DLL_EXPORT
#endif

/* \mainpage TIXI Manual
<b>Contents:</b>
  - \link Glossary Glossary\endlink
  - \link  XPathExamples XPath Examples\endlink
  - \link UsageExamples Usage\endlink
  - <a HREF="modules.html"  class="el">Function Documentation</a>
 */
/*
 @file   tixi.h
 @author Hans-Peter Kersken <Hans-Peter.Kersken@dlr.de>
         Markus Litz <Markus.Litz@dlr.de>,
         Markus Kunde <Markus.Kunde@dlr.de>,
         Arne Bachmann <Arne.Bachmann@dlr.de>,
 @date   Tue Feb 02 14:28:05 2009

 @brief  Definition of  <b>TIXI</b> the <b>TI</b>VA - <b>XML</b> - <b>I</b>nterface
*/
/*! \page Glossary Glossary

 XML-file: file containing an XML-document

 XML-document: abstraction which either refers to the contents of an
 XML-file or to the data structure held in memory.
*/
/*! \page XPathExamples XPath Examples

 The string to describe a path to an element has to be an XPath
 expression but only a very restricted subset of XPath expressions
 are allowed: On retrieving information the specified path has to be
 unique. If an element name appears more than once inside the same
 parent element (see wings/wing in the example) the expression has
 to be contain an index to denominate exactly which element has to
 be used. As long as these restrictions are met any syntactically
 correct XPath expression is allowed. When adding an element the
 complete path has to be specified.

 Example:
 @verbatim
 <?xml version="1.0" encoding="UTF-8"?>
 <plane>
     <name>Junkers JU 52</name>
     <wings>
         <wing position='left'>
             <centerOfGravity referenceSystem='relative'>
                 <x unit="m">30.0</x>
                 <y unit="m">10.0</y>
                 <z unit="m">5.0</z>
             </centerOfGravity>
         </wing>
         <wing position='right'>
             <centerOfGravity referenceSystem='relative'>
                 <x unit="m">30.0</x>
                 <y unit="m">-10.0</y>
                 <z unit="m">5.0</z>
             </centerOfGravity>
         </wing>
     </wings>
     <coordinateOrigin>
         <x unit="m">0.0</x>
         <y unit="m">0.0</y>
         <z unit="m">0.0</z>
     </coordinateOrigin>
 </plane>
@endverbatim

 The path to the x coordinate value of the center of gravity of the
 plane is:

 /plane/centerOfGravity/x

 The path to the x coordinate value of the center of gravity of the
 first wing element is:

 /plane/wings/wing[1]/centerOfGravity/x

 The path to the name can either by expressed as

 /plane/name

 or

 //name

 which would select all elements called name but as long as this element name is
 unique in the document it is valid in the context of TIXI.

 None valid expressions are:

 /plane/wings/wing/centerOfGravity/z is not unique because there are two wing elements in wings.

 /plane/coordinateOrigin specifies an element with children containing text.

 //x is not unique because it would select all elements named x as there are:
  -  /plane/wings/wing[1]/centerOfGravity/x
  -  /plane/wings/wing[2]/centerOfGravity/x
  -  /plane/centerOfGravity/x
*/
/**
 @page UsageExamples Usage

 The following exmaple assume the XML example from section \ref XPathExamples to be stored in a file named ju.xml.

 - Retrieve the position of the first wing.

@verbatim
 #include "tixi.h"

 char* xmlFilename = "ju.xml";
 TixiDocumentHandle handle = -1;
 char * elementPath = "/plane/wings/wing[1]";
 char* attributeName = "position";

 char* attributeValue;

 tixiOpenDocument( xmlFilename, &handle );
 tixiGetTextAttribute( handle, elementPath, attributeName, &attributeValue );
 tixiCloseDocument( handle );
@endverbatim

 - Retrieve x value of the coordinate origin.

@verbatim
 #include "tixi.h"

 char* xmlFilename = "ju.xml";
 TixiDocumentHandle handle = -1;
 char * elementPath = "/plane/coordinateOrigin/x";
 double x = 0.;

 tixiOpenDocument( xmlFilename, &handle );
 tixiGetDoubleElement( handle, elementPath, &x );
 tixiCloseDocument( handle );
@endverbatim
*/
/**
 \page Fortran Notes for Fortran programming

  The Fortran interface is implemented by calls to subroutines.

  It assumes the follwing mapping for the basic types:

  real is real*8 corresponds double

  integer is integer*4 corresponds to int

  character corresponds char

  Character strings are to be passed as variables of type
  character*N. If a string is returned by a subroutine call the
  variable holding the result must be large enough to hold the result.
  Otherwise the result is truncated and the return code is set to
  STRING_TRUNCATED.

  NOTE: In view of these restrictions an implementation using
  character arrays should be considered.

  The return codes returned in the last argument corresponds to their
  position in ::ReturnCode starting with 0. A routine will be supplied
  to directly access the meaning of a return code as a string.

  When the C interface requires to pass a NULL-pointer, e.g. to choose
  the default format "%g" when writing floating point elements, the
  respective argument in the Fortran interfaces is the empty string
  constant "". This is the only way to represent a string of length
  zero. Passing a variable with all characters set to "" will, via the
  interface transformed into an emtpy C-string "\0", which is of
  length 1, and not to a NULL-pointer.

*/


/**
 @internal
 @section ImplementationIssuues Implementation issues

  - Memory allocated by TIXI and associated with a document is released
   when closing the document.

  - TixiDocumentHandle is an integer used as index to access an TIXI
   internal data structure.

 @todo Add for final version:
 @todo - addNamespace
 @todo - check for attribute/elemet name not NULL in arguments passed to respective routines
 @todo - Fortran90, FORTRAN77 interface testing and return code checking

 @todo - when adding elements, created necessary parent elements which do not already exist too.
*/

#ifndef TIXI_H
#define TIXI_H

/**
   Datatype for TixiDocumentHandle.
*/
typedef int TixiDocumentHandle;


/**
  \defgroup Enums Enumerations
*/

/**
    \ingroup Enums
    ReturnCode return code of TIXI-routines.
     Has a typedef to ReturnCode.
*/
enum ReturnCode
{
  SUCCESS,                    /*!< 0: No error occurred                         */

  FAILED,                     /*!< 1: Unspecified error                         */

  INVALID_XML_NAME,           /*!< 2: Non XML standard compliant name specified */

  NOT_WELL_FORMED,            /*!< 3: Document is not well formed               */

  NOT_SCHEMA_COMPLIANT,       /*!< 4: Document is not schema compliant          */

  NOT_DTD_COMPLIANT,          /*!< 5: Document is not DTD compliant             */

  INVALID_HANDLE,             /*!< 6: Document handle is not valid              */

  INVALID_XPATH,              /*!< 7: XPath expression is not valid             */

  ELEMENT_NOT_FOUND,          /*!< 8: Element does not exist in document        */

  INDEX_OUT_OF_RANGE,         /*!< 9: Index supplied as argument is not
                                   inside the admissible range               */

  NO_POINT_FOUND,             /*!< 10: No point element found a given XPath      */

  NOT_AN_ELEMENT,             /*!< 11: XPath expression does not point to an
                                   XML-element node                          */

  ATTRIBUTE_NOT_FOUND,        /*!< 12: Element does not have the attribute       */

  OPEN_FAILED,                /*!< 13: Error on opening the file                 */

  OPEN_SCHEMA_FAILED,         /*!< 14: Error on opening the schema file          */

  OPEN_DTD_FAILED,            /*!< 15: Error on opening the DTD file             */

  CLOSE_FAILED,               /*!< 16: Error on closing the file                 */

  ALREADY_SAVED,              /*!< 17: Trying to modify already saved document   */

  ELEMENT_PATH_NOT_UNIQUE,    /*!< 18: Path expression can not be resolved
                                   unambiguously                             */

  NO_ELEMENT_NAME,            /*!< 19: Element name argument is NULL             */

  NO_CHILDREN,                /*!< 20: Node has no children                      */

  CHILD_NOT_FOUND,            /*!< 21: Named child is not child of element
                                   specified                                 */

  EROROR_CREATE_ROOT_NODE,    /*!< 22: Error when adding root node to new
                                   document                                  */

  DEALLOCATION_FAILED,        /*!< 23: On closing a document the
                                   deallocation of allocated memory fails    */

  NO_NUMBER,                  /*!< 24: No number specified                       */

  NO_ATTRIBUTE_NAME,          /*!< 25: No attribute name specified               */

  STRING_TRUNCATED,           /*!< 26: String variable supplied is to small to
                                   hold the result, Fortran only             */

  NON_MATCHING_NAME,          /*!< 27: Row or column name specified do not
                                   match the names used in the document      */

  NON_MATCHING_SIZE,          /*!< 28: Number of rows or columns specified do
                                   not match the sizes of the matrix in the
                                   document                                  */

  MATRIX_DIMENSION_ERROR,     /*!< 29: if nRows or nColumns or both are
                                   less than 1 */

  COORDINATE_NOT_FOUND,       /*!< 30: missing coordinate inside a point element */

  UNKNOWN_STORAGE_MODE,        /*!< 31: storage mode specified is neither
                                   ROW_WISE nor COLUMN_WISE  */

  UID_NOT_UNIQUE,                /*!< 32: One or more uID's are not unique */

  UID_DONT_EXISTS,                /*!< 33: A given uID's does not exist */

  UID_LINK_BROKEN                /*!< 33: A node the is specified as a Link has no correspoding uid in that data set */

};

typedef enum ReturnCode ReturnCode;

/**

 \ingroup Enums
      Strorage mode of arrays This enum indicates how a matrix is stored
      in an array on input or should be retrieved into an array on
      output. If ROW_WISE is specified the order will be (1,1), (1,2),
      ..., (2,1), (2,2), (2,3)... If COLUMN_WISE is specified the order
      will be (1,1), (2,1), (3,1), .., (1,2),(2,2),(3,2) ...

  Has a typedef to StorageMode.
   */
enum StorageMode
{
  ROW_WISE,                   /*!< row wise order                             */
  COLUMN_WISE                 /*!< column wise                                */
};


typedef enum StorageMode StorageMode;


/**

 \ingroup Enums
      Open mode of a xml file
      This enum indicates how a xml file is opend.
      If OpenMode is OPENMODE_PLAIN, the xml file is open "normal" and just the given
      file is opend. If OpenMode is OPENMODE_RECURSIVE, then all external files
      specified in a <externaldata> node are opend and replaced in the xml tree.

  Has a typedef to OpenMode.
   */
enum OpenMode
{
  OPENMODE_PLAIN,                   /*!< Open just the xml file       */
  OPENMODE_RECURSIVE                /*!< Open with external files     */
};


typedef enum OpenMode OpenMode;

/**

  \ingroup Enums
       Type of a TiXI message

       This is used, if the user wants to override the default message system
       in TiXI.
    */
enum MessageType
{
  MESSAGETYPE_ERROR,                 /*!< The message is an error      */
  MESSAGETYPE_WARNING,               /*!< The message is a warning     */
  MESSAGETYPE_STATUS                 /*!< A status message             */
};


typedef enum MessageType MessageType;

/**
 * TixiPrintMsgFnc:
 * @param[in]  type The message type (error, warning, status)
 * @param[in]  msg The message
 * @param[in]  ... extra arguments
 *
 * Signature of a callback function to handle messages (errors, warnings ...)
 * To be used in conjuction with ::tixiSetPrintMsgFunc.
 */
typedef void (*TixiPrintMsgFnc) (MessageType type, const char *msg);

/**
  @brief Returns the version number of this TIXI version.

  <b>Fortran syntax:</b>
  tixi_get_version( character version )

  @return
      - char* A string with the version number.
*/
DLL_EXPORT char* tixiGetVersion();



/**
  \defgroup FileHandling File Handling Functions
    Function to open, create, close, and save XML-files.
 */
/*@{*/
/**
  @brief Open an XML-file for reading.

  Opens an XML-file specified by xmlFilename for reading and checks
  if it is well formed. To validate the document against a XML-Schema
  or DTD use ::tixiSchemaValidateFromFile, ::tixiSchemaValidateFromString or ::tixiDTDValidate.


  <b>Fortran syntax:</b>

  tixi_open_document( character*n xml_filename, integer handle, integer error )

  @param[in]  xmlFilename name of the XML-file to be opened
  @param[out] handle  handle to the XML-document. This handle is used in
                      calls to other TIXI functions.
  @return

    - SUCCESS if successfully opened the XML-file
    - NOT_WELL_FORMED if opening the XML-file succeeds but test for
                      well-formedness fails
    - OPEN_FAILED if opening of the XML-file failed
 */

DLL_EXPORT ReturnCode tixiOpenDocument (const char *xmlFilename, TixiDocumentHandle * handle);



/**
  @brief Open an XML-file for reading. It acts like tixiOpenDocument.

  If OpenMode is OPENMODE_PLAIN, only the given xml is opened. If OpenMode is
  OPENMODE_RECURSIVE, external xml files will be integrated into the xml tree
  by linking them into the main xml file. The user now has only one big xml file.
  In the path node has be be a valid URI. This uri could adress a relativ or absolut file path,
  or a http url. Example values for the path node are:
    - absolute local directory: "file:///tmp/" or "file:///c:/windws/"
    - relative local directory: "file://relativeDirectory/" or "file://../anotherRelativeDirectory/"
    - remote http ressource: "http://www.someurl.de/"

  Examples for the externaldata node:
  @code{.xml}
  <wings>
      <airfoils>
              <externaldata>
                  <path>file://aDirectory/</path>
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
  @endcode

  <b>Fortran syntax:</b>

  tixi_open_document_recursive( character*n xml_filename, integer handle, integer openmode, integer error )

  @param[in]  xmlFilename name of the XML-file to be opened
  @param[out] handle      handle to the XML-document. This handle is used in
                          calls to other TIXI functions.
  @param[in]  oMode       Enum of the mode to open (OPENMODE_PLAIN / OPENMODE_RECURSIVE).

  @return
    - SUCCESS if successfully opened the XML-file
    - NOT_WELL_FORMED if opening the XML-file succeeds but test for
                      well-formedness fails
    - OPEN_FAILED if opening of the XML-file failed
 */
DLL_EXPORT ReturnCode tixiOpenDocumentRecursive (const char *xmlFilename, TixiDocumentHandle * handle, OpenMode oMode);

/**
  @brief Close an XML-document.

  Closes an XML-document. This routine should be called after the
  processing of an XML-document is completed. After calling this
  routine the handle is invalid and no further processing of the
  document is possible. It must be called after ::tixiSaveDocument.

  <b>Fortran syntax:</b>

  tixi_close_document( integer  handle, integer error )

  @param[in]  handle file handle as returned by ::tixiOpenDocument, ::tixiOpenDocumentRecursive, ::tixiOpenDocumentFromHTTP, ::tixiCreateDocument or ::tixiImportFromString

  @return
    - SUCCESS if successfully closed the XML-file
    - INVALID_HANDLE if  handle not found in list of man
    - CLOSE_FAILED if closing  the XML-file failed
 */

DLL_EXPORT ReturnCode tixiCloseDocument (TixiDocumentHandle handle);

/**
  \defgroup Elements Element handling functions
  Functions to get the content of an element as a string or a number, functions
  to create and manipulate its content, and a funtion to remove an element are described
  in this section.
 */

/*@{*/
/**
  @brief Retrieve text content of an element.

  Returns the text content of the element specified by elementPath in the
  document specified by handle. elementPath must refer to exactly one
  element which has only a text node and zero or more attributes but
  no further children with text nodes. If an error occurs text is set
  to NULL. On successful return the memory used for text is allocated
  internally and must not be released by the user. The deallocation
  is handle when the document referred to by handle is closed.

  <b>Fortran syntax:</b>

  tixi_get_text_element( integer handle, character*n element_path,  character*n text, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP

  @param[in]  elementPath an XPath compliant path to an element in the document
                    specified by handle (see section \ref XPathExamples above).

  @param[out] text text content of the element specified by elementPath

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetTextElement (const TixiDocumentHandle handle,
                                          const char *elementPath, char **text);


/**
  @brief Retrieve integer content of an element.

  Returns the content of the element specified by elementPath in the
  document specified by handle as an integer. elementPath must refer to exactly one
  element which has only a text node and zero or more attributes but
  no further children with text nodes. If an error occurs text is set
  to NULL. On successful return the memory used for text is allocated
  internally and must not be released by the user. The deallocation
  is handle when the document referred to by handle is closed.

  <b>Fortran syntax:</b>

  tixi_get_integer_element( integer  handle, character*n element_path, int* number, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP

  @param[in]  elementPath an XPath compliant path to an element in the document
                    specified by handle (see section \ref XPathExamples above).

  @param[out] number  content of the element specified by elementPath interpreted as an integer number

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetIntegerElement (const TixiDocumentHandle handle, const char *elementPath, int *number);

/**
  @brief Retrieve floating point content of an element.

  Returns the content of the element specified by elementPath in the
  document specified by handle as a floating point
  number. elementPath must refer to exactly one element which has
  only a text node and zero or more attributes but no further
  children with text nodes. If an error occurs number is set to
  NULL.

  <b>Fortran syntax:</b>

  tixi_get_double_element( integer  handle, character*n element_path, real number, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP

  @param[in]  elementPath an XPath compliant path to an element in the document
                    specified by handle (see section \ref XPathExamples above).

  @param[out] number  content of the element specified by elementPath interpreted as a floating point number

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetDoubleElement (const TixiDocumentHandle handle, const char *elementPath, double *number);


/**
  @brief Returns the number of children elements with the same name.

  <b>Fortran syntax:</b>

  tixi_get_named_children_count( integer handle, character*n element_path, character*n child_name, int* count, integer error )

  @param[in]  handle handle as returned by ::tixiCreateDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP
  @param[in]  elementPath the path to an element in the document
                          specified by handle (see section \ref XPathExamples above).
  @param[in]  childName name of children to be counted
  @param[out] count  number of children with name childName.
                     0 is returned if either the element specified by elementPath has no
                     children at all or has no children with name childName.

  @return
    - SUCCESS if a count is computed
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but
                              to a list of elements
    - NO_CHILD_NAME if childName is NULL
 */
DLL_EXPORT ReturnCode tixiGetNamedChildrenCount (const TixiDocumentHandle handle,
                                                 const char *elementPath, const char *childName,
                                                 int *count);

/**
  @brief Returns the name of a child node beneath a given path.

  <b>Fortran syntax:</b>

  tixi_get_child_node_name( integer handle, character*n element_path, int* index, character*n child_name_array, integer error )

  @param[in]  handle handle as returned by ::tixiCreateDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP
  @param[in]  parentElementPath the path to the parent element in the document
                                specified by handle (see section \ref XPathExamples above).
  @param[in]  index number index of the child-element of the given path.
  @param[out] name  String containing the name of the child node. If the node is not a normal node, the name variable will contain:
                     - \#text - in case of a text node
                     - \#comment - in case of a comment node
                     - \#cdata-section - in case of a CDATA section node

  @return
    - SUCCESS if a count is computed
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but
                              to a list of elements
*/
DLL_EXPORT ReturnCode tixiGetChildNodeName (const TixiDocumentHandle handle,
                                            const char *parentElementPath, int index, char **name);

/**
  @brief Returns the number of child elements beneath a given path.

  <b>Fortran syntax:</b>

  tixi_get_number_of_childs( integer handle, character*n element_path, int* nchilds, integer error )

  @param[in]  handle handle as returned by ::tixiCreateDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP
  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).
  @param[out] nChilds Number of child elements beneath the given elementPath.

  @return
    - SUCCESS if a count is computed
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but
                              to a list of elements
*/
DLL_EXPORT ReturnCode tixiGetNumberOfChilds(const TixiDocumentHandle handle, const char *elementPath, int* nChilds);

/*@}*/

/**
  \defgroup Attributes Attribute Handling Functions

  Functions to get the content of an element attribute as a string or a number,
  functions to create and manipulate attributes, and a function to remove attributes
  are described in this section.
 */
/*@{*/

/**
  @brief Retrieves value of an element's attribute as a string.

  Returns the value of an attribute specified by attributeName of the
  element, specified by elementPath, in the document specified by
  handle. On successful return the memory used for value is allocated
  internally and must not be released by the user. The memory is
  deallocated when the document referred to by handle is closed.

  <b>Fortran syntax:</b>

  tixi_get_text_attribute( integer  handle, character*n element_path, character*n attribute_name, character*n text, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument or ::tixiCreateDocument

  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).

  @param[in]  attributeName name of the attribute to be get from the element

  @param[out] text value of the specified attribute as a string

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ATTRIBUTE_NOT_FOUND if the element has no attribute attributeName
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetTextAttribute (const TixiDocumentHandle handle,
                                            const char *elementPath, const char *attributeName,
                                            char **text);

/**
  @brief Retrieves value of an element's attribute as an integer.

  Returns the value of an attribute specified by attributeName of the
  element, specified by elementPath, in the document specified by
  handle. On successful return the memory used for value is allocated
  internally and must not be released by the user. The memory is
  deallocated when the document referred to by handle is closed.

  <b>Fortran syntax:</b>

  tixi_get_integer_attribute( integer  handle, character*n element_path, character*n attribute_name, integer *number, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument or ::tixiCreateDocument

  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).

  @param[in]  attributeName name of the attribute to be added to the element

  @param[out] number  value of the specified attribute as an integer value

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ATTRIBUTE_NOT_FOUND if the element has no attribute attributeName
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetIntegerAttribute (const TixiDocumentHandle handle,
                                               const char *elementPath, const char *attributeName,
                                               int *number);

/**
  @brief Retrieves value of an element's attribute as a floating point number.

  Returns the value of an attribute specified by attributeName of the
  element, specified by elementPath, in the document specified by
  handle. On successful return the memory used for value is allocated
  internally and must not be released by the user. The memory is
  deallocated when the document referred to by handle is closed.

  <b>Fortran syntax:</b>

  tixi_get_double_attribute( integer  handle, character*n element_path, character*n attribute_name, real *number, integer error )

  @param[in]  handle handle as returned by ::tixiOpenDocument or ::tixiCreateDocument

  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).

  @param[in]  attributeName name of the attribute to be added to the element

  @param[out] number value of the specified attribute as a floating point value

  @return
    - SUCCESS if successfully retrieve the text content of a single element
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ATTRIBUTE_NOT_FOUND if the element has no attribute attributeName
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetDoubleAttribute (const TixiDocumentHandle handle,
                                              const char *elementPath, const char *attributeName,
                                              double *number);
                                              
/**
  @brief Returns the number of attributes  of a given node.

  <b>Fortran syntax:</b>

  tixi_get_number_of_attributes( integer handle, character*n element_path, int* nattr, integer error )

  @param[in]  handle handle as returned by ::tixiCreateDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP
  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).
  @param[out] nAttributes Number of attributes of a given node.

  @return
    - SUCCESS if a count is computed
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but
                              to a list of elements
 */
DLL_EXPORT ReturnCode tixiGetNumberOfAttributes(const TixiDocumentHandle handle, const char *elementPath, int* nAttributes);

/**
  @brief Returns the name of an attribute beneath a given path.

  <b>Fortran syntax:</b>

  tixi_get_attribute_name( integer handle, character*n element_path, int* index, character*n attr_name_array, integer error )

  @param[in]  handle handle as returned by ::tixiCreateDocument, ::tixiOpenDocumentRecursive or ::tixiOpenDocumentFromHTTP
  @param[in]  elementPath an XPath compliant path to an element in the document
                          specified by handle (see section \ref XPathExamples above).
  @param[in]  attrIndex number index of the attribute of the given path (counting from 1...tixiGetNumberOfAttributes)
  @param[out] attrName  String containing the attribute name.

  @return
    - SUCCESS if a count is computed
    - INVALID_HANDLE if the handle is not valid, i.e.  does not or no longer exist
    - INVALID_XPATH if elementPath is not a well-formed XPath-expression
    - ELEMENT_NOT_FOUND if elementPath does not point to a node in the XML-document
    - ELEMENT_PATH_NOT_UNIQUE if elementPath resolves not to a single element but
                              to a list of elements
 */

DLL_EXPORT ReturnCode tixiGetAttributeName(const TixiDocumentHandle handle, const char *elementPath, int attrIndex, char** attrName);

/**
  @brief Returns the file path to the document.

  The path is empty, if the document was not opened, but created by ::tixiCreateDocument.

  @param[in]  handle document handle as returned by ::tixiOpenDocument, ::tixiOpenDocumentRecursive, ::tixiOpenDocumentFromHTTP, ::tixiCreateDocument or ::tixiImportFromString
  @param[out] documentPath  Path to the file, opened by ::tixiOpenDocument, ::tixiOpenDocumentRecursive, ::tixiOpenDocumentFromHTTP
                            The path is a null pointer, if the document was created by ::tixiCreateDocument or ::tixiImportFromString

  @return
    - SUCCESS in case of no errors.
    - FAILED if documentPath is a null pointer.
    - INVALID_HANDLE if the document handle is invalid.
 */

DLL_EXPORT ReturnCode tixiGetDocumentPath (TixiDocumentHandle handle, char** documentPath);

/*@}*/


#endif  /* TIXI_H */


#ifdef __cplusplus
}
#endif
