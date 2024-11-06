/*************************
  Copyright 1999, Greg Landrum

	this has the files for doing file opening/choosing on the Macintosh
	
	created by greg Landrum  January 22, 1996
	
*************************/
#include "viewkel.h"
#include "Mac_Fopen.h"
#include <StandardFile.h>

extern FILE *fopen_mac ( short vRefNum , long parID , char * fileName , char * mode,
												char action );

FILE *choose_mac_file(char *file_name,char action)
{
	StandardFileReply reply;
	SFTypeList		type_list;
	short num_types;
	FILE *the_file;
	
	/* start out by using StandardGetFile to get the file name and stuff */
	num_types = -1;
	type_list[0] = 'YhMp';
	StandardGetFile(nil, num_types, type_list, &reply);
	if (!reply.sfGood)
	{
		strcpy(file_name,"User Cancelled");
		return 0;
	}


	/* now open the file */
	strcpy(file_name,(char *)reply.sfFile.name);
	PtoCstr((unsigned char *)file_name);
	the_file = fopen_mac(reply.sfFile.vRefNum,reply.sfFile.parID,
											 file_name,(char *)"r",action);
	
	return the_file;
}


FILE *
fopen_mac ( short vRefNum , long parID , char * fileName , char * mode,
						char action ) 
{
	short oldVol ;
	short aVol ;
	long aDir , aProc ;
	FILE * ret = NULL ;

	/* change to the proper disk and directory */
  if ( GetVol ( NULL , & oldVol ) ) {
     return NULL ;
  }
  if ( GetWDInfo ( oldVol , & aVol , & aDir , & aProc ) ) {
     return NULL  ;
  }
  if ( HSetVol ( NULL , vRefNum , parID ) ) {
     return NULL ;
  }
		if( action == MAC_FOPEN_OPEN_CD || action == MAC_FOPEN_OPEN_NOCD ){
   	 ret = fopen ( fileName , mode ) ;
   	} else{
   	 ret = 0;
   	}
		
		/* this doesn't seem to work the way I want it to... */
    if( action == MAC_FOPEN_OPEN_NOCD || action == MAC_FOPEN_NOOPEN_NOCD ){
   	  if ( HSetVol ( NULL, aVol , aDir ) ) {
    	    /* an error we can't currently handle */
    	}
    	if ( SetVol ( NULL, oldVol ) ) {
        /* an error we can't currently handle */
    	}
    } 
    return ret ;
}


