/*
* Copyright (C) 20015 German Aerospace Center (DLR/SC)
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
* Altered by Sven Werchner in order to create a minimal version for ICON-ART
* Note: When XML-functionality for ICON-ART is expanding, 
*       this TIXI-Version probably needs to be (re-)expanded as well
*/

#include "uidHelper.h"

extern void printMsg(MessageType type, const char* message, ...);

int uid_clearUIDList(TixiDocument *document)
{
  TixiUIDListEntry *current = document->uidListHead;
  while (current) {
    TixiUIDListEntry *next = (TixiUIDListEntry *) current->next;
    if (current->uIDName) {
      xmlFree(current->uIDName);
      current->uIDName = NULL;
    }
    free(current);
    current = next;
  }
  document->uidListHead = NULL;
  return SUCCESS;
}