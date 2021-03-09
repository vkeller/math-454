MODULE DlgSaveScreen
	USE globals
    !USE i f w i n
    USE ContWrap
    USE ResWrap
    USE CmnDlgs
    USE Printfunc

    PRIVATE
    PUBLIC ScreenToBitmapFile
    SAVE

    
	CHARACTER(LEN=40)				   :: fname
	CHARACTER(LEN=256)				   :: fullpath
	
	!   bitmap handle has module scope so it can
	!   be accessed by the proc function
	INTEGER(HANDLE)                    :: hbitmap, hdib


CONTAINS


    SUBROUTINE ScreenToBitmapFile
		USE charfunc
		USE filesubs
		
		! need to replace CVF interfaces
		!USE GDI32, Ignore_GetDIBits => GetDIBits, Ignore_cp => CreatePalette
		USE GDI32, Ignore_cp => CreatePalette
        
		IMPLICIT NONE

		INTEGER(HANDLE)                     :: hdc, hmemdc, holdbitmap, hpal, ihandl
		INTEGER								:: ncf, ncp
		INTEGER								:: rval
		INTEGER								:: ncolors, lp
		INTEGER								:: infosize, width, height
											
		TYPE(T_BITMAP)						:: bm
		TYPE(T_BITMAPINFOHEADER)			:: bi
		TYPE(T_BITMAPFILEHEADER)			:: bfh
		TYPE(T_RGBQUAD)						:: rgbq
		TYPE(T_PALETTEENTRY),DIMENSION(256) :: pe

		!	this interface is now included in CVF6.6c, 
		!	so a custom interface is no longer needed
		!INTERFACE
		!	FUNCTION GetDIBits (arg1, arg2, arg3, arg4, arg5, arg6, arg7)
		!		USE ifwinTY
		!		integer(SINT) :: GetDIBits ! int
		!		!DEC$ ATTRIBUTES DEFAULT, STDCALL, DECORATE, ALIAS:'GetDIBits' :: GetDIBits
		!		integer(HANDLE) arg1 ! HDC arg1
		!		integer(HANDLE) arg2 ! HBITMAP arg2
		!		integer(UINT) arg3 ! UINT arg3
		!		integer(UINT) arg4 ! UINT arg4
		!		integer(LPVOID) arg5 ! LPVOID arg5
		!		integer(LPVOID) arg6 ! LPBITMAPINFO arg6
		!		integer(UINT) arg7 ! UINT arg7
		!	END FUNCTION
		!END INTERFACE

		!	custom interface to CreatePalette to pass the arg as an integer LOC()
		!	passed by value instead of GDI32's version where the argument is a TYPE
		!	passed by reference
		INTERFACE 
			FUNCTION CreatePalette (arg1)
				USE ifwinTY
				integer(HANDLE) :: CreatePalette ! HPALETTE
				!DEC$ ATTRIBUTES DEFAULT, STDCALL, DECORATE, ALIAS:'CreatePalette' :: CreatePalette
				integer(LPVOID) arg1 ! LPLOGPALETTE arg1
			END FUNCTION
		END INTERFACE

		! Each row in the BMP is padded to next higher 4 byte boundary.
		INTEGER					:: Pad, sx
		Pad(sx) = ISHFT(ISHFT((sx)+3,-2),2)


		!     TYPE T_BITMAP
		!     SEQUENCE
		!       integer(LONG)	 bmType			! knowns  LONG 
		!       integer(LONG)	 bmWidth		! knowns  LONG 
	!       integer(LONG)	 bmHeight		! knowns  LONG 
		!       integer(LONG)	 bmWidthBytes	! knowns  LONG 
		!       integer(WORD)	 bmPlanes		! knowns  WORD 
		!       integer(WORD)	 bmBitsPixel	! knowns  WORD 
		!       integer(LPVOID)  bmBits			! knowns  LPVOID 
		!     END TYPE

		!     TYPE T_RGBQUAD
		!     SEQUENCE
		!       integer(BYTE)	 rgbBlue		! knowns  BYTE 
		!       integer(BYTE)	 rgbGreen		! knowns  BYTE 
		!       integer(BYTE)	 rgbRed			! knowns  BYTE 
		!       integer(BYTE)	 rgbReserved	! knowns  BYTE 
		!     END TYPE

		!	  TYPE T_BITMAPINFOHEADER
		!	  SEQUENCE
		!	    integer(DWORD)	 biSize			! knowns  DWORD 
		!	    integer(LONG)	 biWidth		! knowns  LONG 
		!	    integer(LONG)	 biHeight		! knowns  LONG 
		!	    integer(WORD)	 biPlanes		! knowns  WORD 
		!	    integer(WORD)	 biBitCount		! knowns  WORD 
		!	    integer(DWORD)	 biCompression	! knowns  DWORD 
		!	    integer(DWORD)	 biSizeImage	! knowns  DWORD 
		!	    integer(LONG)	 biXPelsPerMeter ! knowns  LONG 
		!	    integer(LONG)	 biYPelsPerMeter ! knowns  LONG 
		!	    integer(DWORD)	 biClrUsed		! knowns  DWORD 
		!	    integer(DWORD)	 biClrImportant ! knowns  DWORD 
		!	  END TYPE

		!     TYPE T_BITMAPFILEHEADER
		!     SEQUENCE
		!       integer(WORD)	 bfType			! knowns  WORD 
		!       integer(DWORD)	 bfSize			! knowns  DWORD 
		!       integer(WORD)	 bfReserved1	! knowns  WORD 
		!       integer(WORD)	 bfReserved2	! knowns  WORD 
		!       integer(DWORD)	 bfOffBits		! knowns  DWORD 
		!     END TYPE

		!     TYPE T_LOGPALETTE
		!     SEQUENCE
		!       integer(WORD) palVersion ! knowns  WORD 
		!       integer(WORD) palNumEntries ! knowns  WORD 
		!       TYPE (T_PALETTEENTRY) palPalEntry(1) ! typedefs  PALETTEENTRY 
		!     END TYPE


	!======= STEP 1: PREPARE THE BITMAP==========================================
	! this is done BEFORE the savetofile dialog, so that the screen capture is
	! immediate, and not screwed up by any time delays in clearing the savetofile
	! dialog

		!	abstract the entire screen (null GetDC arg) as a bitmap
		hdc     = GetDC                  (NULL)
		width   = GetDeviceCaps          (hdc, HORZRES)
		height  = GetDeviceCaps          (hdc, VERTRES)	
		hbitmap = CreateCompatibleBitmap (hdc, width, height)

		!	memory device context
		hmemdc     = CreateCompatibleDC (hdc)
		holdbitmap = SelectObject       (hmemdc, hbitmap)
		
		!	copy the image data into the memory device context
		rval = BitBlt (hmemdc, 0,0,		&
					   width, height,	&
					   hdc, 0,0,		&
					   SRCCOPY)

		!	instead of using the LOGPALETTE type, simply block out
		!	some memory and accomplish the same thing, only this allows
		!	the paletteentry() array to be associated with a varying 
		!	number of colors
		IF (IAND(GetDeviceCaps (hdc, RASTERCAPS), RC_PALETTE) > 0) THEN
			lp = MALLOC (8 + 256*SIZEOF(pe(1)))
			rval = #300										! palVersion
			CALL CopyMemory (LOC(lp), LOC(rval), 4)
			rval = GetSystemPaletteEntries (hdc, 0, 255, pe(1))
			CALL CopyMemory (LOC(lp)+4, LOC(rval), 4)		! palNumEntries
			CALL CopyMemory (LOC(lp)+8, LOC(pe(1)), 256*SIZEOF(pe(1)))
			hpal = CreatePalette (lp)
			CALL FREE (lp)
		ELSE
			hpal = GetStockObject (DEFAULT_PALETTE)
		END IF
		rval = ReleaseDC (NULL, hdc)

		!	get the bitmap structure members
		rval = GetObject (hbitmap, SIZEOF(bm), LOC(bm))

		!	fill in the bitmap information structure
		bi%biSize			= SIZEOF(bi)
		bi%biWidth			= Pad(bm%bmWidth)
		bi%biHeight			= bm%bmHeight
		bi%biPlanes			= 1
		bi%biBitCount		= bm%bmPlanes * bm%bmBitsPixel
		bi%biCompression	= BI_RGB
		bi%biSizeImage		= 0
		bi%biXPelsPerMeter	= 0
		bi%biYPelsPerMeter	= 0
		bi%biClrUsed		= 0
		bi%biClrImportant	= 0

		!	number of colors
		SELECT CASE (bi%biBitCount)
		CASE (1)
			ncolors = 2
		CASE (4)
			ncolors = 16
		CASE (8)
			ncolors = 256
		CASE DEFAULT
			ncolors = 0
		END SELECT

		!	size of infoheader and color table
		infosize = bi%biSize + ncolors * SIZEOF(rgbq)

		!	create a device context for the DIB
		hdc  = GetDC (NULL)
		hpal = MSFWIN$SelectPalette (hdc, hpal, FALSE)
		rval = RealizePalette (hdc)

		!	allocate memory for the infoheader and color table
		hdib = GlobalAlloc (GMEM_FIXED, infosize)
		IF (hdib == 0) THEN
			rval = MSFWIN$SelectPalette (hdc, hpal, FALSE)
			rval = DeleteObject (hpal)
			rval = ReleaseDC (NULL, hdc)
			RETURN
		END IF
		CALL CopyMemory (hdib, LOC(bi), bi%biSize)

		!	get SizeImage from the device driver
		rval = GetDIBits (hmemdc,			&	! source device context
						  hbitmap,			&	! bitmap handle
						  0,				&	! first scan line
						  bi%biHeight,		&	! number of lines to copy
						  0,				&	! null--> fills in bitmapinfo only
						  hdib,				&	! addr of bitmapInfo structure
						  DIB_RGB_COLORS)		! use RGB colors		   
		
		CALL CopyMemory (LOC(bi), hdib, bi%biSize)
		IF (bi%biSizeImage == 0) THEN
			bi%biSizeIMage = bi%biHeight *	&
							 (IAND((bi%biWidth*bi%biBitCount + 31), .NOT.(31))/8)
		END IF

		!	enlarge the buffer to hold the pixel data
		CALL CopyMemory (hdib, LOC(bi), bi%biSize)
		hdib = GlobalReAlloc (hdib, infosize + bi%biSizeImage, GMEM_MOVEABLE)

		!	get the entire DIB
		rval = GetDIBits (hmemdc,			&	! source device context
						  hbitmap,			&	! bitmap handle
						  0,				&	! first scan line
						  bi%biHeight,		&	! number of lines to copy
						  hdib + infosize,	&	! memory offset to start addr. for pixel data
						  hdib,				&	! addr of bitmapInfo structure
						  DIB_RGB_COLORS)		! use RGB colors		   

		rval = MSFWIN$SelectPalette (hdc, hpal, FALSE)
		rval = ReleaseDC (NULL, hdc)



	!======= STEP 2: SHOW A DIALOG TO GET THE FILENAME & PATH=======================

		CALL concat (RootPath, BmpDir, fullpath)
		fname = 'screen'

		IF (ShowModalDialog(IDD_SAVESCREEN, LOC(ScreenSaverProc)) == IDOK) THEN

			!	prepare the target filename from the provided strings
			ncf = INDEX(fname, '.')
			IF (ncf > 0) fname(ncf:) = ''
			ncf = chcnt (fname, 20)
			ncp = chcnt (fullpath, 200)
			
			IF (ncf > 0 .AND. ncp > 0) THEN
				fullpath = fullpath(1:ncp)//'\'//fname(1:ncf)//'.bmp'//CHAR(0)

				!	create the file
				ihandl = open_the_file (fullpath, 'W')
				IF (ihandl > 0) THEN

					!	bitmap file header
					bfh%bfType		= MakeWord (INT1(ICHAR('B')), INT1(ICHAR('M')))
					bfh%bfSize		= GlobalSize(hdib) + SIZEOF(bfh)
					bfh%bfReserved1 = 0
					bfh%bfReserved2 = 0
					bfh%bfOffBits   = SIZEOF(bfh) + infosize
					CALL rw_file ('W', ihandl, SIZEOF(bfh), LOC(bfh))
					
					!	bitmap info header + colormap + pixel data
					CALL rw_file ('W', ihandl, infosize + bi%biSizeImage, hdib)
					
					CALL close_file (ihandl)
				END IF
			END IF
		END IF

		!	release system resources
		rval = DeleteObject (hpal)
		rval = DeleteObject (hbitmap)
		rval = DeleteDC     (hmemdc)
		rval = GlobalFree   (hdib)

    END SUBROUTINE ScreenToBitmapFile



    INTEGER FUNCTION ScreenSaverProc (hwnd, msg, wParam, lParam) RESULT(res)
    !DEC$ IF DEFINED(_X86_)
    !DEC$ ATTRIBUTES STDCALL, ALIAS : '_ScreenSaverProc@16' :: ScreenSaverProc
    !DEC$ ELSE
    !DEC$ ATTRIBUTES STDCALL, ALIAS : 'ScreenSaverProc' :: ScreenSaverProc
    !DEC$ ENDIF
        IMPLICIT NONE
        INTEGER(HANDLE),  INTENT(IN)   :: hwnd
        INTEGER(UINT),    INTENT(IN)   :: msg
        INTEGER(fWPARAM), INTENT(IN)   :: wParam
        INTEGER(fLPARAM), INTENT(IN)   :: lParam

        INTEGER                        :: rval
        INTEGER(HANDLE)                :: hwndControl
        INTEGER                        :: controlId
        INTEGER                        :: code, print_orientation
        
        INTEGER(HANDLE)                :: hdc

        res = DefaultDialogProc (hwnd, msg, wParam, lParam)

        SELECT CASE (msg)

        CASE (WM_INITDIALOG)
		    CALL EditBoxSetText (hwnd, IDC_SS_PATH, fullpath)        
		    CALL EditBoxSetText (hwnd, IDC_SS_FILE, fname)        
            res = 1

		CASE (WM_PAINT)
			CALL ktlogo_bctr (hwnd)
			res = 1

        CASE (WM_COMMAND)
            code		= HIWORD(wParam)
            controlId   = LOWORD(wParam)
            hwndControl = lParam

            SELECT CASE (controlId)

            CASE (IDOK)
				CALL EditBoxReadText (hwnd, IDC_SS_FILE, fname, 20)
				CALL EditBoxReadText (hwnd, IDC_SS_PATH, fullpath, 200)
				res = 1

            CASE (IDC_SS_BROWSE)
                IF (code == BN_CLICKED) THEN
					IF (CmnDlgChooseFolder (hwnd, fullpath))    &
						CALL EditBoxSetText (hwnd, IDC_SS_PATH, fullpath)
                END IF
                res = 1

			CASE (IDC_SS_PRINT)
				IF (code == BN_CLICKED) THEN
					IF (printer_setup (hwnd, hdc, print_orientation)) THEN
						rval = SetCursor (LoadCursor(NULL, IDC_WAIT))
				        rval = EnableWindow	(hwnd, FALSE)
				        rval = StartPage	(hdc)
			            rval = print_bitmap (hdc, 0, hbitmap,   &	! bitmap handle
								             10, 10,	        &	! origin coords
								             1., top_left)    	    ! sizefactor & origin location
				        rval = EndPage		(hdc)
				        rval = EndDoc		(hdc)
				        rval = EnableWindow	(hwnd, TRUE)
						rval = SetCursor (LoadCursor(NULL, IDC_ARROW))
						rval = SetFocus (hwnd)
					END IF
				END IF
				res = 1
            
            END SELECT

        CASE (WM_DESTROY)

        END SELECT
        
    END FUNCTION ScreenSaverProc



END MODULE DlgSaveScreen
