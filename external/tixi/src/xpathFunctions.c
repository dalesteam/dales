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

#include "xpathFunctions.h"


extern void printMsg(MessageType type, const char* message, ...);

xmlXPathObjectPtr XPathEvaluateExpression(xmlDocPtr doc, const char* xPathExpression)
{
  xmlXPathContextPtr xpathContext;
  xmlXPathObjectPtr xpathObject;

  /* Create xpath evaluation context */
  xpathContext = xmlXPathNewContext(doc);
  if (!xpathContext) {
    printMsg(MESSAGETYPE_ERROR, "Error: Unable to create new XPath context.\n");
    return NULL;
  }

  /* Evaluate Expression */
  xpathObject = xmlXPathEvalExpression((xmlChar*) xPathExpression, xpathContext);
  if (!(xpathObject)) {
    printMsg(MESSAGETYPE_ERROR, "Error: Invalid XPath expression \"%s\"\n", xPathExpression);
    xmlXPathFreeContext(xpathContext);
    return NULL;
  }

  if (xmlXPathNodeSetIsEmpty(xpathObject->nodesetval)) {
    xmlXPathFreeContext(xpathContext);
    xmlXPathFreeObject(xpathObject);
    return NULL;
  }

  xmlXPathFreeContext(xpathContext);
  return xpathObject;
}