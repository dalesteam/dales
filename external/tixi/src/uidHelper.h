/*
* Copyright (C) 2015 German Aerospace Center (DLR/SC)
*
* Created: 2010-08-13 Markus Litz <Markus.Litz@dlr.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*   http://www.apache.org/licenses/LICENSE-2.0
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

/**
 * @file   uidHelper.h
 * @author Markus Litz <Markus.Litz@dlr.de>
 * @date   Tue Jul 25 12:06:33 2010
 *
 * @brief  Header file for helper routines with UIDs.
 *
 */
#include "tixi.h"
#include "tixiData.h"

/**
 * Cleans up memory and removed the uid-list from the document.
 */
int uid_clearUIDList(TixiDocument *document);

