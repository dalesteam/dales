! #############################################################################
! # Copyright (C) 2015 German Aerospace Center (DLR/SC)
! #
! # Created: 2015-09-16 Melven Roehrig-Zoellner <Melven.Roehrig-Zoellner@DLR.de>
! #
! # Licensed under the Apache License, Version 2.0 (the "License");
! # you may not use this file except in compliance with the License.
! # You may obtain a copy of the License at
! #
! #     http://www.apache.org/licenses/LICENSE-2.0
! #
! # Unless required by applicable law or agreed to in writing, software
! # distributed under the License is distributed on an "AS IS" BASIS,
! # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! # See the License for the specific language governing permissions and
! # limitations under the License.
! #############################################################################
! 
! # This file is automatically created from tixi.h on 2016-01-13.
! # If you experience any bugs please contact the authors
! 


!
! Altered by Sven Werchner (KIT) in order to create a minimal version for ICON-ART
! Note: When XML-functionality for ICON-ART is expanding, 
!       this TIXI-Version probably needs to be (re-)expanded as well
!


module tixi
  use, intrinsic :: iso_c_binding, only: C_CHAR
  implicit none
  public


! enum StorageMode
enum, bind(C)
  enumerator :: ROW_WISE = 0
  enumerator :: COLUMN_WISE = 1
end enum

! enum ReturnCode
enum, bind(C)
  enumerator :: SUCCESS = 0
  enumerator :: FAILED = 1
  enumerator :: INVALID_XML_NAME = 2
  enumerator :: NOT_WELL_FORMED = 3
  enumerator :: NOT_SCHEMA_COMPLIANT = 4
  enumerator :: NOT_DTD_COMPLIANT = 5
  enumerator :: INVALID_HANDLE = 6
  enumerator :: INVALID_XPATH = 7
  enumerator :: ELEMENT_NOT_FOUND = 8
  enumerator :: INDEX_OUT_OF_RANGE = 9
  enumerator :: NO_POINT_FOUND = 10
  enumerator :: NOT_AN_ELEMENT = 11
  enumerator :: ATTRIBUTE_NOT_FOUND = 12
  enumerator :: OPEN_FAILED = 13
  enumerator :: OPEN_SCHEMA_FAILED = 14
  enumerator :: OPEN_DTD_FAILED = 15
  enumerator :: CLOSE_FAILED = 16
  enumerator :: ALREADY_SAVED = 17
  enumerator :: ELEMENT_PATH_NOT_UNIQUE = 18
  enumerator :: NO_ELEMENT_NAME = 19
  enumerator :: NO_CHILDREN = 20
  enumerator :: CHILD_NOT_FOUND = 21
  enumerator :: EROROR_CREATE_ROOT_NODE = 22
  enumerator :: DEALLOCATION_FAILED = 23
  enumerator :: NO_NUMBER = 24
  enumerator :: NO_ATTRIBUTE_NAME = 25
  enumerator :: STRING_TRUNCATED = 26
  enumerator :: NON_MATCHING_NAME = 27
  enumerator :: NON_MATCHING_SIZE = 28
  enumerator :: MATRIX_DIMENSION_ERROR = 29
  enumerator :: COORDINATE_NOT_FOUND = 30
  enumerator :: UNKNOWN_STORAGE_MODE = 31
  enumerator :: UID_NOT_UNIQUE = 32
  enumerator :: UID_DONT_EXISTS = 33
  enumerator :: UID_LINK_BROKEN = 34
end enum

! enum MessageType
enum, bind(C)
  enumerator :: MESSAGETYPE_ERROR = 0
  enumerator :: MESSAGETYPE_WARNING = 1
  enumerator :: MESSAGETYPE_STATUS = 2
end enum

! enum OpenMode
enum, bind(C)
  enumerator :: OPENMODE_PLAIN = 0
  enumerator :: OPENMODE_RECURSIVE = 1
end enum




abstract interface
    subroutine TixiPrintMsgFnc(messageType, msg) bind(C)
        use, intrinsic :: iso_c_binding
        integer(kind=C_INT), value :: messageType
        character(kind=C_CHAR), intent(in) :: msg(*)
    end subroutine
end interface




interface
! char* tixiGetVersion();
  function tixiGetVersion_c() &
      result(ret) &
      bind(C,name='tixiGetVersion')
    use, intrinsic :: iso_c_binding
    type(C_PTR) :: ret
  end function tixiGetVersion_c


! ReturnCode tixiOpenDocument (const char *xmlFilename, TixiDocumentHandle * handle);
  function tixiOpenDocument_c(xmlFilename, &
        handle) &
      result(ret) &
      bind(C,name='tixiOpenDocument')
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR), intent(in) :: xmlFilename(*)
    integer(kind=C_INT), intent(out) :: handle
    integer(kind=C_INT) :: ret
  end function tixiOpenDocument_c
  
! ReturnCode tixiOpenDocumentRecursive (const char *xmlFilename, TixiDocumentHandle * handle, OpenMode oMode);
  function tixiOpenDocumentRecursive_c(xmlFilename, &
        handle, &
        oMode) &
      result(ret) &
      bind(C,name='tixiOpenDocumentRecursive')
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR), intent(in) :: xmlFilename(*)
    integer(kind=C_INT), intent(out) :: handle
    integer(kind=C_INT), value :: oMode
    integer(kind=C_INT) :: ret
  end function tixiOpenDocumentRecursive_c


! ReturnCode tixiCloseDocument (TixiDocumentHandle handle);
  function tixiCloseDocument(handle) &
      result(ret) &
      bind(C,name='tixiCloseDocument')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    integer(kind=C_INT) :: ret
  end function tixiCloseDocument

! ReturnCode tixiGetTextElement (const TixiDocumentHandle handle, const char *elementPath, char **text);
  function tixiGetTextElement_c(handle, &
        elementPath, &
        text) &
      result(ret) &
      bind(C,name='tixiGetTextElement')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    type(C_PTR), intent(inout) :: text
    integer(kind=C_INT) :: ret
  end function tixiGetTextElement_c


! ReturnCode tixiGetIntegerElement (const TixiDocumentHandle handle, const char *elementPath, int *number);
  function tixiGetIntegerElement_c(handle, &
        elementPath, &
        number) &
      result(ret) &
      bind(C,name='tixiGetIntegerElement')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    integer(kind=C_INT), intent(out) :: number
    integer(kind=C_INT) :: ret
  end function tixiGetIntegerElement_c


! ReturnCode tixiGetDoubleElement (const TixiDocumentHandle handle, const char *elementPath, double *number);
  function tixiGetDoubleElement_c(handle, &
        elementPath, &
        number) &
      result(ret) &
      bind(C,name='tixiGetDoubleElement')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    real(kind=C_DOUBLE), intent(out) :: number
    integer(kind=C_INT) :: ret
  end function tixiGetDoubleElement_c

! ReturnCode tixiGetNamedChildrenCount (const TixiDocumentHandle handle, const char *elementPath, const char *childName, int *count);
  function tixiGetNamedChildrenCount_c(handle, &
        elementPath, &
        childName, &
        count) &
      result(ret) &
      bind(C,name='tixiGetNamedChildrenCount')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    character(kind=C_CHAR), intent(in) :: childName(*)
    integer(kind=C_INT), intent(out) :: count
    integer(kind=C_INT) :: ret
  end function tixiGetNamedChildrenCount_c


! ReturnCode tixiGetChildNodeName (const TixiDocumentHandle handle, const char *parentElementPath, int index, char **name);
  function tixiGetChildNodeName_c(handle, &
        parentElementPath, &
        index, &
        name) &
      result(ret) &
      bind(C,name='tixiGetChildNodeName')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: parentElementPath(*)
    integer(kind=C_INT), value :: index
    type(C_PTR), intent(inout) :: name
    integer(kind=C_INT) :: ret
  end function tixiGetChildNodeName_c


! ReturnCode tixiGetNumberOfChilds(const TixiDocumentHandle handle, const char *elementPath, int* nChilds);
  function tixiGetNumberOfChilds_c(handle, &
        elementPath, &
        nChilds) &
      result(ret) &
      bind(C,name='tixiGetNumberOfChilds')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    integer(kind=C_INT), intent(out) :: nChilds
    integer(kind=C_INT) :: ret
  end function tixiGetNumberOfChilds_c


! ReturnCode tixiGetTextAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, char **text);
  function tixiGetTextAttribute_c(handle, &
        elementPath, &
        attributeName, &
        text) &
      result(ret) &
      bind(C,name='tixiGetTextAttribute')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    character(kind=C_CHAR), intent(in) :: attributeName(*)
    type(C_PTR), intent(inout) :: text
    integer(kind=C_INT) :: ret
  end function tixiGetTextAttribute_c


! ReturnCode tixiGetIntegerAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, int *number);
  function tixiGetIntegerAttribute_c(handle, &
        elementPath, &
        attributeName, &
        number) &
      result(ret) &
      bind(C,name='tixiGetIntegerAttribute')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    character(kind=C_CHAR), intent(in) :: attributeName(*)
    integer(kind=C_INT), intent(out) :: number
    integer(kind=C_INT) :: ret
  end function tixiGetIntegerAttribute_c

! ReturnCode tixiGetDoubleAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, double *number);
  function tixiGetDoubleAttribute_c(handle, &
        elementPath, &
        attributeName, &
        number) &
      result(ret) &
      bind(C,name='tixiGetDoubleAttribute')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    character(kind=C_CHAR), intent(in) :: attributeName(*)
    real(kind=C_DOUBLE), intent(out) :: number
    integer(kind=C_INT) :: ret
  end function tixiGetDoubleAttribute_c

! ReturnCode tixiGetNumberOfAttributes (const TixiDocumentHandle handle, const char *elementPath, int nattr);
  function tixiGetNumberOfAttributes_c(handle, &
        elementPath, &
        nattr) &
      result(ret) &
      bind(C,name='tixiGetNumberOfAttributes')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    integer(kind=C_INT), intent(out) :: nattr
    integer(kind=C_INT) :: ret
  end function tixiGetNumberOfAttributes_c

! ReturnCode tixiGetAttributeName (const TixiDocumentHandle handle, const char *elementPath, int *attrIndex, const char *attrName);
  function tixiGetAttributeName_c(handle, &
        elementPath, &
        attrIndex, &
        text_c) &
      result(ret) &
      bind(C,name='tixiGetAttributeName')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    character(kind=C_CHAR), intent(in) :: elementPath(*)
    integer(kind=C_INT), value :: attrIndex
    type(C_PTR), intent(inout) :: text_c
    integer(kind=C_INT) :: ret
  end function tixiGetAttributeName_c

! ReturnCode tixiGetDocumentPath (TixiDocumentHandle handle, char** documentPath);
  function tixiGetDocumentPath_c(handle, &
        text_c) &
      result(ret) &
      bind(C,name='tixiGetDocumentPath')
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), value :: handle
    type(C_PTR), intent(inout) :: text_c
    integer(kind=C_INT) :: ret
  end function tixiGetDocumentPath_c


end interface


  private :: tixiGetVersion_c
  private :: tixiOpenDocument_c
  private :: tixiOpenDocumentRecursive_c
  private :: tixiGetTextElement_c
  private :: tixiGetIntegerElement_c
  private :: tixiGetDoubleElement_c
  private :: tixiGetNamedChildrenCount_c
  private :: tixiGetChildNodeName_c
  private :: tixiGetNumberOfChilds_c
  private :: tixiGetTextAttribute_c
  private :: tixiGetIntegerAttribute_c
  private :: tixiGetDoubleAttribute_c
  private :: tixiGetNumberOfAttributes_c
  private :: tixiGetAttributeName_c
  private :: tixiGetDocumentPath_c

  private :: f_c_strarrayptr
  private :: c_f_strarrayptr
  private :: c_f_stringptr


! we need a wrapper type to create arrays of C-strings as there is no pointer of a pointer in Fortran
type CStringPtr
  character(kind=C_CHAR), pointer :: str(:)
end type


contains
  
  subroutine f_c_strarrayptr(strarray_f,tmp_str,tmp_ptr,strarray_c)
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR,len=*), intent(in) :: strarray_f(:)
    character(kind=C_CHAR), target, allocatable, intent(inout) :: tmp_str(:,:)
    type(C_PTR), target, allocatable, intent(inout) :: tmp_ptr(:)
    type(C_PTR), intent(inout) :: strarray_c
    integer :: i, j
  
    allocate(tmp_str(1+len(strarray_f),size(strarray_f)))
    allocate(tmp_ptr(size(strarray_f)))
    do i = 1, size(strarray_f), 1
      do j = 1, len(strarray_f), 1
        tmp_str(j,i) = strarray_f(i)(j:j)
      end do
      tmp_str(len(strarray_f)+1,i) = C_NULL_CHAR
      tmp_ptr(i) = c_loc( tmp_str(1,i) )
    end do
  
    strarray_c = c_loc(tmp_ptr)
  end subroutine
  
  ! WARNING this function is written such that it
  !         circumvents triggering old GCC gfortran bugs
  subroutine c_f_strarrayptr(n,strarray_c,strarray_f)
    use, intrinsic :: iso_c_binding
    integer, intent(in) :: n
    type(C_PTR), intent(in) :: strarray_c(n)
    type(CStringPtr), target, intent(inout) :: strarray_f(n)
    type(CStringPtr), pointer :: tmp => null()
    integer :: i
  
    do i = 1, n, 1
      tmp => strarray_f(i)
      call c_f_stringptr(strarray_c(i),tmp%str)
    end do
  end subroutine
  
  subroutine c_f_stringptr(str_c, str_f)
    use, intrinsic :: iso_c_binding
    type(C_PTR), intent(in) :: str_c
    character(kind=C_CHAR), pointer, intent(inout) :: str_f(:)
    integer :: i
  
    if( .not. c_associated(str_c) ) return
    i = 1
    call c_f_pointer(str_c,str_f,(/i/))
    do while(str_f(i) .ne. C_NULL_CHAR)
      i = i + 1
      call c_f_pointer(str_c,str_f,(/i/))
    end do
    call c_f_pointer(str_c,str_f,(/i-1/))
  
  end subroutine c_f_stringptr




  
  subroutine tixiGetVersion(str)
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=C_CHAR), pointer, intent(out) :: str(:)
    type(C_PTR) :: str_c
  
    str_c = tixiGetVersion_c()
    call c_f_stringptr(str_c,str)
  
  end subroutine
  
  ! commented out for gfortran <= 4.3 compatibility
  !function tixiGetPrintMsgFunc() result(fnc)
  !  use, intrinsic :: iso_c_binding
  !  implicit none
  !  procedure(TixiPrintMsgFnc), pointer :: fnc
  !  type(C_FUNPTR) :: fnc_c
  !  type(C_PTR) :: fnc_c_
  !
  !  fnc_c_ = tixiGetPrintMsgFunc_c()
  !  ! ugly conversion from C_PTR to C_FUNPTR
  !  fnc_c = transfer(fnc_c_, fnc_c)
  !  call c_f_procpointer(fnc_c,fnc)
  !end function
  !
  !
  !function tixiSetPrintMsgFunc(fnc) result(ret)
  !  use, intrinsic :: iso_c_binding
  !  implicit none
  !  procedure(TixiPrintMsgFnc), pointer, intent(in) :: fnc
  !  integer(kind=C_INT) :: ret
  !  type(C_FUNPTR) :: fnc_c
  !  type(C_PTR) :: fnc_c_
  !
  !  fnc_c = c_funloc(fnc)
  !  ! ugly conversion from C_FUNPTR to C_PTR
  !  fnc_c_ = transfer(fnc_c, fnc_c_)
  !  ret = tixiSetPrintMsgFunc_c(fnc_c_)
  !end function




! ReturnCode tixiOpenDocument (const char *xmlFilename, TixiDocumentHandle * handle);
  function tixiOpenDocument(xmlFilename, &
        handle) &
      result(ret)
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR,len=*), intent(in) :: xmlFilename
    integer(kind=C_INT), intent(out) :: handle
    integer(kind=C_INT) :: ret

    ret = tixiOpenDocument_c(xmlFilename // C_NULL_CHAR, &
        handle)

  end function tixiOpenDocument


! ReturnCode tixiOpenDocumentRecursive (const char *xmlFilename, TixiDocumentHandle * handle, OpenMode oMode);
  function tixiOpenDocumentRecursive(xmlFilename, &
        handle, &
        oMode) &
      result(ret)
    use, intrinsic :: iso_c_binding
    character(kind=C_CHAR,len=*), intent(in) :: xmlFilename
    integer(kind=C_INT), intent(out) :: handle
    integer(kind=C_INT), intent(in) :: oMode
    integer(kind=C_INT) :: ret

    ret = tixiOpenDocumentRecursive_c(xmlFilename // C_NULL_CHAR, &
        handle, &
        oMode)

  end function tixiOpenDocumentRecursive


! ReturnCode tixiGetTextElement (const TixiDocumentHandle handle, const char *elementPath, char **text);
  function tixiGetTextElement(handle, &
        elementPath, &
        text) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    character(kind=C_CHAR), pointer :: text(:)
    integer(kind=C_INT) :: ret
    type(C_PTR) :: text_c = C_NULL_PTR


    ret = tixiGetTextElement_c(handle, &
        elementPath // C_NULL_CHAR, &
        text_c)

    call c_f_stringptr(text_c, text)

  end function tixiGetTextElement


! ReturnCode tixiGetIntegerElement (const TixiDocumentHandle handle, const char *elementPath, int *number);
  function tixiGetIntegerElement(handle, &
        elementPath, &
        number) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    integer(kind=C_INT), intent(out) :: number
    integer(kind=C_INT) :: ret

    ret = tixiGetIntegerElement_c(handle, &
        elementPath // C_NULL_CHAR, &
        number)

  end function tixiGetIntegerElement


! ReturnCode tixiGetDoubleElement (const TixiDocumentHandle handle, const char *elementPath, double *number);
  function tixiGetDoubleElement(handle, &
        elementPath, &
        number) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    real(kind=C_DOUBLE), intent(out) :: number
    integer(kind=C_INT) :: ret

    ret = tixiGetDoubleElement_c(handle, &
        elementPath // C_NULL_CHAR, &
        number)

  end function tixiGetDoubleElement


! ReturnCode tixiGetNamedChildrenCount (const TixiDocumentHandle handle, const char *elementPath, const char *childName, int *count);
  function tixiGetNamedChildrenCount(handle, &
        elementPath, &
        childName, &
        count) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    character(kind=C_CHAR,len=*), intent(in) :: childName
    integer(kind=C_INT), intent(out) :: count
    integer(kind=C_INT) :: ret

    ret = tixiGetNamedChildrenCount_c(handle, &
        elementPath // C_NULL_CHAR, &
        childName // C_NULL_CHAR, &
        count)

  end function tixiGetNamedChildrenCount


! ReturnCode tixiGetChildNodeName (const TixiDocumentHandle handle, const char *parentElementPath, int index, char **name);
  function tixiGetChildNodeName(handle, &
        parentElementPath, &
        index, &
        name) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: parentElementPath
    integer(kind=C_INT), intent(in) :: index
    character(kind=C_CHAR), pointer :: name(:)
    integer(kind=C_INT) :: ret
    type(C_PTR) :: name_c = C_NULL_PTR


    ret = tixiGetChildNodeName_c(handle, &
        parentElementPath // C_NULL_CHAR, &
        index, &
        name_c)

    call c_f_stringptr(name_c, name)

  end function tixiGetChildNodeName


! ReturnCode tixiGetNumberOfChilds(const TixiDocumentHandle handle, const char *elementPath, int* nChilds);
  function tixiGetNumberOfChilds(handle, &
        elementPath, &
        nChilds) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    integer(kind=C_INT), intent(out) :: nChilds
    integer(kind=C_INT) :: ret

    ret = tixiGetNumberOfChilds_c(handle, &
        elementPath // C_NULL_CHAR, &
        nChilds)

  end function tixiGetNumberOfChilds


! ReturnCode tixiGetTextAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, char **text);
  function tixiGetTextAttribute(handle, &
        elementPath, &
        attributeName, &
        text) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    character(kind=C_CHAR,len=*), intent(in) :: attributeName
    character(kind=C_CHAR), pointer :: text(:)
    integer(kind=C_INT) :: ret
    type(C_PTR) :: text_c = C_NULL_PTR


    ret = tixiGetTextAttribute_c(handle, &
        elementPath // C_NULL_CHAR, &
        attributeName // C_NULL_CHAR, &
        text_c)

    call c_f_stringptr(text_c, text)

  end function tixiGetTextAttribute


! ReturnCode tixiGetIntegerAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, int *number);
  function tixiGetIntegerAttribute(handle, &
        elementPath, &
        attributeName, &
        number) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    character(kind=C_CHAR,len=*), intent(in) :: attributeName
    integer(kind=C_INT), intent(out) :: number
    integer(kind=C_INT) :: ret

    ret = tixiGetIntegerAttribute_c(handle, &
        elementPath // C_NULL_CHAR, &
        attributeName // C_NULL_CHAR, &
        number)

  end function tixiGetIntegerAttribute

! ReturnCode tixiGetDoubleAttribute (const TixiDocumentHandle handle, const char *elementPath, const char *attributeName, double *number);
  function tixiGetDoubleAttribute(handle, &
        elementPath, &
        attributeName, &
        number) &
      result(ret)
    use, intrinsic :: iso_c_binding
    integer(kind=C_INT), intent(in) :: handle
    character(kind=C_CHAR,len=*), intent(in) :: elementPath
    character(kind=C_CHAR,len=*), intent(in) :: attributeName
    real(kind=C_DOUBLE), intent(out) :: number
    integer(kind=C_INT) :: ret

    ret = tixiGetDoubleAttribute_c(handle, &
        elementPath // C_NULL_CHAR, &
        attributeName // C_NULL_CHAR, &
        number)

  end function tixiGetDoubleAttribute

! ReturnCode tixiGetNumberOfAttributes(const TixiDocumentHandle handle, const char *elementPath, int* nAttributes);
  function tixiGetNumberOfAttributes(handle, &
         elementPath,  &
         nattr) &
    result(ret)
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), intent(in) :: handle
  character(kind=C_CHAR,len=*), intent(in) :: elementPath
  integer(kind=C_INT), intent(out) :: nattr
  integer(kind=C_INT) :: ret

  ret = tixiGetNumberOfAttributes_c(handle, &
       elementPath // C_NULL_CHAR, &
       nattr)
  end function tixiGetNumberOfAttributes


!ReturnCode tixiGetAttributeName(const TixiDocumentHandle handle, const char *elementPath, int attrIndex, char** attrName);
  function tixiGetAttributeName(handle, &
         elementPath,  &
         attr_idx, &
         attrName) &
    result(ret)
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), intent(in) :: handle
  character(kind=C_CHAR,len=*), intent(in) :: elementPath
  integer(kind=C_INT), intent(in) :: attr_idx
  character(kind=C_CHAR), pointer :: attrName(:)
  integer(kind=C_INT) :: ret
  type(C_PTR) :: text_c = C_NULL_PTR

  ret = tixiGetAttributeName_c(handle, &
       elementPath // C_NULL_CHAR, &
       attr_idx,  &
       text_c)

    call c_f_stringptr(text_c, attrName)
  end function tixiGetAttributeName


!ReturnCode tixiGetDocumentPath (TixiDocumentHandle handle, char** documentPath);
  function tixiGetDocumentPath(handle, &
         docPath) &
    result(ret)
  use, intrinsic :: iso_c_binding
  integer(kind=C_INT), intent(in) :: handle
  character(kind=C_CHAR), pointer :: docPath(:)
  integer(kind=C_INT) :: ret
  type(C_PTR) :: text_c = C_NULL_PTR

  ret = tixiGetDocumentPath_c(handle, &
       text_c)

    call c_f_stringptr(text_c, docPath)
  end function tixiGetDocumentPath



end module
