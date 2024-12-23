/*******************************************************
*      Copyright (C) 1996, 1999 Greg Landrum
*
*  This file is part of yaehmop.
*
*   This is free software.
* 
*  Permission is granted to modify, or otherwise fold, spindle, and mutilate this
*    code provided all copyright notices are left intact.
*
*  This code may be distributed to your heart's content, in whatever form,
*    provided no fee is charged for the distribution, all copyright notices are
*    left intact, and the source is distributed (without fee) along with any
*    binaries to anyone who requests it.
*
*  There are, of course, no warranties at all on this program.
*
********************************************************************/

#define CLOSE_TO_ZERO 1e-4

/********

  this has got the stuff for dealing with contour plots....

*********/  

/***
  Recent Edit History

   26.09.98 gL:
     various modifications to remove warning when compiled with -Wall
     under gcc (yeah yeah... it's anal)
   22.01.99 gL:
     changes to contouring code so that it works better with non-FCO
     plots (specifically so that it can read in stuff generated by
     the dump_grids functionality from MO_conts.c).
***/
  

#include "viewkel.h"



/****************************************************************************
 *
 *                   Procedure contour_data_set
 *
 * Arguments: 
 *            
 * Returns: none
 *
 * Action:  Does the contouring necessary for a data set.
 *    
 *  
 ****************************************************************************/
void contour_data_set(int num_args,char **cont_plot_ptr)
{
  char instring[80];
  contour_plot_type *cont_plot;
  gnuplot_contour_type *cont_ptr;
  int i,j;
  float *max_x,*max_y,*min_x,*min_y;
  float *old_max_x,*old_max_y,*old_min_x,*old_min_y;
  float last_val;
  int old_cont_kind,num_cont_vals;

  cont_plot = (contour_plot_type *)cont_plot_ptr[0];

  max_x = &(cont_plot->max_x);
  max_y = &(cont_plot->max_y);
  min_x = &(cont_plot->min_x);
  min_y = &(cont_plot->min_y);
  old_max_x = &(cont_plot->old_max_x);
  old_max_y = &(cont_plot->old_max_y);
  old_min_x = &(cont_plot->old_min_x);
  old_min_y = &(cont_plot->old_min_y);
  
  /* check to see if we need to re-determine the location of tic marks */
  if( *old_max_x != *max_x ||
     *old_max_y != *max_y ||
     *old_min_x != *min_x ||
     *old_min_y != *min_y ){


    preprocess_cont_plot_data(cont_plot);
    
    *old_max_x = *max_x;
    *old_max_y = *max_y;
    *old_min_x = *min_x;
    *old_min_y = *min_y;
  }
  


  /**************

    contour the data.  

  ***************/
  
  /* if there is currently some contour data, blow it out first. */
  for(i=0;i<cont_plot->num_curves;i++){
    while(cont_plot->contours[i]){
      cont_ptr = cont_plot->contours[i];
      cont_plot->contours[i] = cont_plot->contours[i]->next;
      D_FREE(cont_ptr->coords);
      D_FREE(cont_ptr);
    }
    cont_plot->contours[i] = 0;
  }

  /********

    we may need some input from the user about the values of the contours

  ********/
  if( cont_plot->levels_kind == LEVELS_INCREMENTAL ){
    display("Look in the xterm");
    if( !cont_plot->levels_list ){
      cont_plot->levels_list = (float *)D_CALLOC(2,sizeof(float));
      if( !(cont_plot->levels_list) ) fatal("Can't allocate levels list");

      /* set some defaults */
      cont_plot->levels_list[0] = 0.01;
      cont_plot->levels_list[1] = 0.01;
    }
    printf("\nEnter the lower limit on the contours [%lf]: ",cont_plot->levels_list[0]);
    fflush(stdin);
    fgets(instring,80,stdin);
    if(instring[0] != '\n' && instring[0] != 0){
      cont_plot->levels_list[0] = (float)atof(instring);
    }
    printf("\nEnter the increment for the contours [%lf]: ",cont_plot->levels_list[1]);
    fflush(stdin);
    fgets(instring,80,stdin);
    if(instring[0] != '\n' && instring[0] != 0){
      cont_plot->levels_list[1] = (float)atof(instring);
    }

    
  }

  /* now contour each curve in the data set */
  fprintf(stderr,"Contouring: \n");
  for(i=0;i<cont_plot->num_curves;i++){
    fprintf(stderr,"%d \n",i);
    cont_plot->contours[i] = contour_data(cont_plot->num_x,i,cont_plot->data,
					  cont_plot->num_levels,
					  cont_plot->num_approx_pts,
					  cont_plot->interp_kind,
					  cont_plot->order,
					  cont_plot->levels_kind,
					  cont_plot->levels_list);

    /********

      use the contours from the first data set to
      determine all the contour levels

    *********/
    if( !i && cont_plot->type == CONT_FCO ){
      cont_ptr = cont_plot->contours[i];
      last_val = -10;
      num_cont_vals = 0;
      while(cont_ptr){
	if( cont_ptr->coords[0].z != last_val ){
	  num_cont_vals++;
	}
	last_val = cont_ptr->coords[0].z;
	cont_ptr = cont_ptr->next;
      }
      cont_plot->levels_list = (double *)D_CALLOC(num_cont_vals,sizeof(double));
      if(!cont_plot->levels_list) fatal("Can't get levels_list");

      fprintf(stderr,"Contour Values:\n");
      j=0;
      cont_ptr = cont_plot->contours[i];
      last_val = -10;
      while(cont_ptr){
	if( cont_ptr->coords[0].z != last_val ){
	  cont_plot->levels_list[j++] = cont_ptr->coords[0].z;
	  fprintf(stderr,"%6.3lf\n",cont_ptr->coords[0].z);
	}
	last_val = cont_ptr->coords[0].z;
	cont_ptr = cont_ptr->next;
      }
      old_cont_kind = cont_plot->levels_kind;
      cont_plot->levels_kind = LEVELS_DISCRETE;
      cont_plot->num_levels = num_cont_vals;
    }
      
      
    if( cont_plot->type != CONT_FCO ){
      fprintf(stderr,"Values:\n");
      cont_ptr = cont_plot->contours[i];
      last_val = -10;
      while(cont_ptr){
	if( cont_ptr->coords[0].z != last_val )
	  fprintf(stderr,"%6.3lf\n",cont_ptr->coords[0].z);
	last_val = cont_ptr->coords[0].z;
      
	cont_ptr = cont_ptr->next;
      }
    }
  }

  if( cont_plot->type == CONT_FCO ){
    cont_plot->levels_kind = old_cont_kind;
  }
}



/****************************************************************************
 *
 *                   Procedure cont_find_tic_sep
 *
 * Arguments: cont_plot: pointer to contour_plot_type
 *            
 * Returns: none
 *
 * Action:  Determines a "reasonable" tic mark separation for the data set
 *    stored in 'graph.
 *  
 ****************************************************************************/
void cont_find_tic_sep(contour_plot_type *cont_plot)
{
  float range,norm,tics,tic_sep,log_range;
  float xmax,xmin,ymax,ymin;

  /* first do the x tics */
  range = fabs(cont_plot->min_x-cont_plot->max_x);
  
  log_range = log10(range);

  norm = pow(10.0,log_range-(float)((log_range >= 0.0 ) ? (int)log_range :
				    ((int)log_range-1)));
  if (norm <= 2)
    tics = 0.2;
  else if (norm <= 5)
    tics = 0.5;
  else tics = 1.0;	
  tic_sep = tics * pow(10.0,(float)((log_range >= 0.0 ) ? (int)log_range : ((int)log_range-1)));

  cont_plot->tic_sep_x = tic_sep;

  /* figure out the max and min values spanned by the tic marks */
  xmin = tic_sep * floor(cont_plot->min_x/tic_sep);
  xmax = tic_sep * ceil(cont_plot->max_x/tic_sep);
  cont_plot->num_tics_x = 1 + (xmax - xmin) / tic_sep;

  /* if the first tic mark is outside the range the user specified, skip it! */
  if( (int)(10000.0*xmin) < (int)(10000.0*cont_plot->min_x) ){
    cont_plot->tic_start_x = xmin + tic_sep;
    xmin += tic_sep;
    cont_plot->num_tics_x--;
  }
  else{
    cont_plot->tic_start_x = xmin;
  }
  /* same deal with the last tic mark */
  if( (int)(10000.0*xmax) > (int)(10000.0*cont_plot->max_x) ){
    cont_plot->num_tics_x--;
    xmax -= tic_sep;
  }
  if( (int)(10000.0*xmax) > (int)(10000.0*cont_plot->max_x) ){
    cont_plot->num_tics_x--;
    xmax -= tic_sep;
  } 

  /******

    repeat the process for the y data

  *******/
  range = fabs(cont_plot->min_y-cont_plot->max_y);
  log_range = log10(range);
  norm = pow(10.0,log_range-(float)((log_range >= 0.0 ) ? (int)log_range :
				    ((int)log_range-1)));
  if (norm <= 2)
    tics = 0.2;
  else if (norm <= 5)
    tics = 0.5;
  else tics = 1.0;	
  tic_sep = tics * pow(10.0,(float)((log_range >= 0.0 ) ? (int)log_range :
				    ((int)log_range-1)));
  cont_plot->tic_sep_y = tic_sep;
  ymin = tic_sep * floor(cont_plot->min_y/tic_sep);
  ymax = tic_sep * ceil(cont_plot->max_y/tic_sep);
  cont_plot->num_tics_y = 1 + (ymax - ymin) / tic_sep;
  if( ymin < cont_plot->min_y ){
    cont_plot->tic_start_y = ymin + tic_sep;
    cont_plot->num_tics_y--;
    ymin += tic_sep;
  }
  else{
    cont_plot->tic_start_y = ymin;
  }
  if( ymax > cont_plot->max_y ){
    cont_plot->num_tics_y--;
  }
/*  if( ymax > cont_plot->max_y ){
    cont_plot->num_tics_y--;
  }
*/
}

/****************************************************************************
 *
 *                   Procedure read_3D_data
 *
 * Arguments: infile: pointer to type FILE
 *               cont_plot: a pointer to contour_plot_type
 *            
 * Returns: none
 *
 * Action:  Reads the data out of 'infile.
 *       This reads until it hits a line beginning with a # mark.
 *
 ****************************************************************************/
void read_3D_data(FILE *infile,contour_plot_type *cont_plot)
{
  char instring[MAX_STR_LEN];
  char tempstring[MAX_STR_LEN];
  int num_curves;
  int xcnt,ycnt;
  int i;
  int max_p,num_p;
  float xval,yval;
  float min_z,max_z,tempval;
  iso_curve_type *temp_iso_curve;
  point_type *temp_point;
  
  if( !cont_plot )FATAL_BUG("No cont_plot passed to read_3D_data.");

  if( skipcomments(infile,instring) < 0){
    error("That file is empty!");
    display("Too Bad....");
    return;
  }

  /*******

    make sure that we've read past any header information

    write some code to process this stuff!

  ********/
  upcase(instring);
  while(instring[0] == '#' && !strstr(instring,"BEGIN_DATA")){
    if( strstr(instring,"NUM_CURVES") ){
      sscanf(instring,"%s %d",tempstring,&cont_plot->num_curves);
      num_curves = cont_plot->num_curves;
    }
    else if( strstr(instring,"FCO_DATA") ){
      cont_plot->type = CONT_FCO;
    }
    else if( strstr(instring,"MO_DATA") ){
      cont_plot->type = CONT_MO;
    }
    else if( strstr(instring,"MIN_X") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->min_x);
    }
    else if( strstr(instring,"MAX_X") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->max_x);
    }
    else if( strstr(instring,"NUM_X") ){
      sscanf(instring,"%s %d",tempstring,&cont_plot->raw_num_x);
    }
    else if( strstr(instring,"MIN_Y") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->min_y);
    }
    else if( strstr(instring,"MAX_Y") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->max_y);
    }
    else if( strstr(instring,"NUM_Y") ){
      sscanf(instring,"%s %d",tempstring,&cont_plot->raw_num_y);
    }
    else if( strstr(instring,"STEP_X") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->step_x);
    }
    else if( strstr(instring,"STEP_Y") ){
      sscanf(instring,"%s %lf",tempstring,&cont_plot->step_y);
    }
    
    if(skipcomments(infile,instring)<0){
      error("That file is empty!");
      display("Too Bad....");
      return;
    }
    upcase(instring);
  }
  if(num_curves < 1){
    error("That data file has less than one data set in it... check it please.");
    display("Yikes!");
    return;
  }

  /* get space to keep track of which curves should be displayed */
  cont_plot->curves_to_display = (char *)D_CALLOC(num_curves,sizeof(char));
  if( !cont_plot->curves_to_display )
    fatal("Can't get space for curves_to_display.");

  /* get space to store linestyles */
  cont_plot->styles = (char *)D_CALLOC(num_curves,sizeof(char));
  if( !cont_plot->styles )fatal("Can't get space for curve styles.");

  /* initialize the styles */
  for( i=0;i<num_curves; i++){
    cont_plot->styles[i] = i;
  }

  /* default to showing just the first curve */
  cont_plot->curves_to_display[0] = 1;

  cont_plot->do_x_tics = 1;
  cont_plot->do_y_tics = 1;

  /*********

    allocate and read in the grid of data

    to make life easier in the contouring functions, we're gonna
    store this data differently than the way it is done in the
    other graphing functions.  For the 3D data, we'll use
    gnuplot style isolines.

    each isoline contains all the curve data, but each
    curve is stored sequentially.   i.e. it's all the
    points for curve 1, followed by all the points for curve 2, etc.

  **********/
  cont_plot->raw_data = 0;
  for(i=0;i<cont_plot->raw_num_x;i++){
    /* the iso curves are stored as a linked list, allocate that now */
    temp_iso_curve = (iso_curve_type *)D_CALLOC(1,sizeof(iso_curve_type));
    if( !temp_iso_curve ){
      error("Can't get space for an iso_curve");
      display("Darn!");
      return;
    }
    temp_iso_curve->next = cont_plot->raw_data;
    cont_plot->raw_data = temp_iso_curve;
    temp_iso_curve->points = (point_type *)D_CALLOC(cont_plot->raw_num_y*num_curves,
						  sizeof(point_type));
    if( !temp_iso_curve->points ){
      error("Can't get space for iso_curve points");
      display("Darn!");
      return;
    }
  }

  min_z = 1e10;
  max_z = -1e10;


  temp_iso_curve = cont_plot->raw_data;
  max_p = cont_plot->raw_num_y*num_curves*cont_plot->raw_num_x;
  num_p = 0;
  skipcomments(infile,instring);
  xval = cont_plot->min_x;
  for(xcnt=0;xcnt<cont_plot->raw_num_x;xcnt++){
    yval = cont_plot->min_y;
    for(ycnt=0;ycnt<cont_plot->raw_num_y;ycnt++){

      /* use strtok to read out space delimited numbers */
      temp_point = &(temp_iso_curve->points[ycnt]);
      sscanf((const char *)strtok(instring," "),"%lf",
	     &(temp_point->z));
      temp_point->x = xval;
      temp_point->y = yval;
      if(temp_point->z > max_z) max_z = temp_point->z;
      if(temp_point->z < min_z) min_z = temp_point->z;

      for(i=1;i<num_curves;i++){
	temp_point = &(temp_iso_curve->points[i*cont_plot->raw_num_y+ycnt]);
	sscanf((const char *)strtok(0," "),"%lf",
	       &(temp_point->z));
	temp_point->x = xval;
	temp_point->y = yval;
	if(temp_point->z > max_z) max_z = temp_point->z;
	if(temp_point->z < min_z) min_z = temp_point->z;

      }

      if( skipcomments(infile,instring) < 0 )
	fatal("EOF hit while reading grid data");
      yval += cont_plot->step_y;
      num_p++;
      if( num_p > max_p ){
	fatal("Too many points were read.  Halting read.");
	xcnt = cont_plot->raw_num_x;
	ycnt = cont_plot->raw_num_y;
      }
    }
    xval += cont_plot->step_x;
    temp_iso_curve->p_max = temp_iso_curve->p_count = cont_plot->raw_num_y;
    temp_iso_curve = temp_iso_curve->next;
  }
  cont_plot->num_p = num_p;
  cont_plot->max_z = max_z;
  cont_plot->min_z = min_z;
  if( cont_plot->max_x < cont_plot->min_x ){
    tempval = cont_plot->min_x;
    cont_plot->min_x = cont_plot->max_x;
    cont_plot->max_x = tempval;
  }
  if( cont_plot->max_y < cont_plot->min_y ){
    tempval = cont_plot->min_y;
    cont_plot->min_y = cont_plot->max_y;
    cont_plot->max_y = tempval;
  }

  /*******

    okay, if we're doing an FCO plot, read out the total and fragment
    DOSs now

    we are guaranteed the same number of points and same
     ranges in each direction....

  *******/
  if( cont_plot->type == CONT_FCO ){
    cont_plot->max_DOS = 0;
    while(instring[0] != '#' && !strstr(instring,"BEGIN_DOS") ){
      if( skipcomments(infile,instring) < 0 )
	fatal("Can't read DOS data for the FCO plot");
      upcase(instring);
    }
    /* get memory for the DOSs */
    cont_plot->raw_data2D =
      (point_type2D *)D_CALLOC((unsigned)cont_plot->raw_num_x*(num_curves+1),
			     sizeof(point_type2D));
    if( !cont_plot->raw_data2D ) fatal("Can't get space for cont_plot 2D data.");

    xval = cont_plot->min_x;
    for(xcnt=0;xcnt<cont_plot->raw_num_x;xcnt++){
      if(skipcomments(infile,instring)<0) fatal("Can't read 2D data");

      sscanf((const char *)strtok(instring," "),"%lf",
	     &(cont_plot->raw_data2D[xcnt*(num_curves+1)].y));
      for(i=1;i<=num_curves;i++){
	sscanf((const char *)strtok(0," "),"%lf",
	       &(cont_plot->raw_data2D[xcnt*(num_curves+1)+i].y));
	if( cont_plot->raw_data2D[xcnt*(num_curves+1)+i].y >
	   cont_plot->max_DOS ){
	  cont_plot->max_DOS = cont_plot->raw_data2D[xcnt*(num_curves+1)+i].y;
	}
      }
      for(i=0;i<=num_curves;i++){
	cont_plot->raw_data2D[xcnt*(num_curves+1)+i].x = xval;
      }
      
      xval += cont_plot->step_x;
    }
  }

  /* get space for the contour array */
  cont_plot->contours = (gnuplot_contour_type **)
    D_CALLOC(cont_plot->num_curves,sizeof(gnuplot_contour_type *));
  if( !cont_plot->contours ) fatal("Can't get space for contour array");
      
  /* that's it */
}


/****************************************************************************
 *
 *                   Procedure preprocess_cont_plot_data
 *
 * Arguments: cont_plot: a pointer to cont_plot_type
 *            
 * Returns: none
 *
 * Action: This does any preprocessing which is needed on the cont_plot data.
 *  
 ****************************************************************************/
void preprocess_cont_plot_data(contour_plot_type *cont_plot)
{
  float xscale,yscale;
  int i,j,jtab;
  int num_along_y;
  int num_curves;
  iso_curve_type *raw_isocurve,*isocurve,*temp_iso_curve;
  point_type *point;
  double xval,yval;
  int num_tot,num_frag;
  int begin_y;
  float cont_xscale,cont_yscale;
  
  /* find the tic mark spacing */
  cont_find_tic_sep(cont_plot);

  num_curves = cont_plot->num_curves;

  /* scale all the points to fit in a DEF_CONT_PLOT_X x DEF_CONT_PLOT_Y box */
  xscale = DEF_CONT_PLOT_X/(cont_plot->max_x - cont_plot->min_x);
  yscale = DEF_CONT_PLOT_Y/(cont_plot->max_y - cont_plot->min_y);
  cont_yscale = (cont_plot->old_max_y - cont_plot->old_min_y) /
    (cont_plot->max_y - cont_plot->min_y);
  cont_xscale = (cont_plot->old_max_x - cont_plot->old_min_x) /
    (cont_plot->max_x - cont_plot->min_x);

  /* if we currently have some data, free it up now */
  while(cont_plot->data){
    temp_iso_curve = cont_plot->data;
    cont_plot->data = cont_plot->data->next;
    D_FREE(temp_iso_curve->points);
    D_FREE(temp_iso_curve);
  }

  cont_plot->data = 0;

  /**********

    cull the 3D points, ditching any which are outside the range
    and scaling those that are in range.

  **********/

  /* start out by finding the first valid raw curve */
  raw_isocurve = cont_plot->raw_data;
  while(raw_isocurve->points[0].x < cont_plot->min_x )

    raw_isocurve = raw_isocurve->next;

  /* now loop through the valid curves */
  cont_plot->num_x = 0;
  num_along_y= -1;
  while(raw_isocurve && raw_isocurve->points[0].x <= cont_plot->max_x &&
	raw_isocurve->points[0].x >= cont_plot->min_x){
    xval = raw_isocurve->points[0].x;
    cont_plot->num_x++;

    /* allocate space for this isosurface */
    isocurve = (iso_curve_type *)D_CALLOC(1,sizeof(iso_curve_type));
    if( !isocurve ){
      error("Can't get space for an iso_curve in preprocessing");
      display("Darn!");
      return;
    }
    isocurve->next = cont_plot->data;
    cont_plot->data = isocurve;
    
    /******

      we don't really know how much memory will be needed along the y direction,
      so we have to figure it out the first time through, from then on things will
      be fine.  

      Note that this is making the assumption that all of the iso lines are
      of the same length. 

    *******/
    if( num_along_y == -1 ){
      /* count the points we'll be using along the y direction */
      num_along_y = 0;
      for(i=0;i<cont_plot->raw_num_y;i++){
	yval = raw_isocurve->points[i].y; 
	if( (yval <= cont_plot->max_y && yval >= cont_plot->min_y) ){
	  if( num_along_y == 0 ) begin_y = i;
	  num_along_y++;
	}
      }
      cont_plot->num_y = num_along_y;
    }

    /* allocate memory for the isocurve now */
    isocurve->points = (point_type *)D_CALLOC(num_along_y*num_curves,sizeof(point_type));
    if( !isocurve->points ) {
      error("Can't get space for iso_curve points in preprocessing");
      fprintf(stderr,"\t I was trying to get: %d\n",num_along_y*num_curves);
      display("Darn!");
      return;
    }

    /* copy in the scaled data */
    for(j=0;j<num_curves;j++){
      jtab = j*cont_plot->num_y;
      for(i=0;i<num_along_y;i++){    
	yval = raw_isocurve->points[j*(cont_plot->raw_num_y)+i+begin_y].y;
	point = &isocurve->points[jtab+i];
	point->x = xval*xscale;
	point->y = yval*yscale;
	point->z = raw_isocurve->points[j*cont_plot->raw_num_y+i+begin_y].z;
      }
    }
    isocurve->p_max = isocurve->p_count = num_along_y;
    raw_isocurve = raw_isocurve->next;
  }
  
  /* adjust the tic mark spacing to fit the new scaling */
  cont_plot->tic_sep_x *= xscale;
  cont_plot->tic_sep_y *= yscale;  
  cont_plot->tic_start_x *= xscale;
  cont_plot->tic_start_y *= yscale;


  if( cont_plot->type == CONT_FCO ){
    /********************

      for FCO plots, the total DOS and fragment DOSs are drawn
      along the Y and X axes respectively.  raw_data2D has
      the unscaled data with the total DOS in the first slot,
      then the fragment DOSs.  The X coordinate provides energy,
      the Y coordinate provides the value.  We'll show the whole
      DOS of each one, in a box with the same X (fragment) or
      Y (total) window and scale as the contour plot.  The DOS value
      will be scaled to fit in a box CONT_PLOT_PROPS_SCALE the
      size of the contour plot window.

    *********************/

    if( !cont_plot->data2D){
      cont_plot->data2D = (point_type2D *)
	D_CALLOC(cont_plot->raw_num_x*(cont_plot->num_curves+1),
	       sizeof(point_type2D));
      if(!cont_plot->data2D ) fatal("Can't allocate cont_plot->data2D");
    }
    num_tot = 0;
    num_frag = 0;
    for( i=0; i<cont_plot->raw_num_x; i++ ){

      /* do the total DOS first */
      xscale = DEF_CONT_PLOT_Y / 
	(cont_plot->max_y - cont_plot->min_y);
      yscale = DEF_CONT_PLOT_X*CONT_PLOT_PROPS_SCALE/
	(cont_plot->max_DOS);
      

#if 0
      if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)].x >
	 cont_plot->max_y ){
      }
      else if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)].x <
	      cont_plot->min_y ){
      }
      else{
#endif
	cont_plot->data2D[i*(cont_plot->num_curves+1)].x =
	  cont_plot->raw_data2D[i*(cont_plot->num_curves+1)].x * xscale;
	if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)].y
	   <= cont_plot->max_DOS ){
	  cont_plot->data2D[i*(cont_plot->num_curves+1)].y =
	    cont_plot->raw_data2D[i*(cont_plot->num_curves+1)].y *
	      yscale;
	}else{
	  cont_plot->data2D[i*(cont_plot->num_curves+1)].y =
	    cont_plot->max_DOS * yscale;
	}
	num_tot++;
#if 0
      }
#endif      
      /* now do the fragment DOSs */
      xscale = DEF_CONT_PLOT_X / 
	(cont_plot->max_x - cont_plot->min_x);
      yscale = DEF_CONT_PLOT_Y*CONT_PLOT_PROPS_SCALE/
	(cont_plot->max_DOS);

#if 0
      if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)+1].x >
	 cont_plot->max_x ){
      }
      else if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)+1].x <
	      cont_plot->min_x ){
      }
      else{
#endif
	for(j=1;j<=cont_plot->num_curves;j++){
	  cont_plot->data2D[i*(cont_plot->num_curves+1)+j].x =
	    cont_plot->raw_data2D[i*(cont_plot->num_curves+1)+j].x * xscale;
	  if( cont_plot->raw_data2D[i*(cont_plot->num_curves+1)+j].y
	     <= cont_plot->max_DOS ){
	    cont_plot->data2D[i*(cont_plot->num_curves+1)+j].y =
	      cont_plot->raw_data2D[i*(cont_plot->num_curves+1)+j].y *
		yscale;
	  }else{
	    cont_plot->data2D[i*(cont_plot->num_curves+1)+j].y =
	      cont_plot->max_DOS * yscale;
	  }
	}
#if 0
	num_frag++;
      }
#endif
    }
  }
#if 0
  /* if there is currently some contour data, blow it out first. */
  for(i=0;i<cont_plot->num_curves;i++){
    while(cont_plot->contours[i]){
      cont_ptr = cont_plot->contours[i];
      cont_plot->contours[i] = cont_plot->contours[i]->next;
      D_FREE(cont_ptr);
    }
    cont_plot->contours[i] = 0;
  }


  /* now contour each curve in the data set */
  fprintf(stderr,"Contouring: \n");
  for(i=0;i<cont_plot->num_curves;i++){
    fprintf(stderr,"%d \n",i);
    cont_plot->contours[i] = contour_data(cont_plot->num_x,i,cont_plot->data,
					  cont_plot->num_levels,
					  cont_plot->num_approx_pts,
					  cont_plot->interp_kind,
					  cont_plot->order,
					  cont_plot->levels_kind,
					  cont_plot->levels_list);
  }
#endif

}
    

/****************************************************************************
 *
 *                   Procedure new_cont_plot
 *
 * Arguments: filename: pointer to type char (optional)
 *            
 * Returns: none
 *
 * Action: does everything to get space for and read in a new cont_plot
 *
 ****************************************************************************/
void new_cont_plot(char *filename)
{
  char file_name[80];
  char *theinline;
  char failed;
  FILE *infile;
  
  failed = 0;
  /* set up a new object to hold the cont_plot */
  makenewobject();
  whichobj = head->obj;

  /* now build the cont_plot primitive */
  whichobj->prim = (prim_type *)D_CALLOC(1,sizeof(prim_type));
  if( !whichobj->prim )fatal("Can't get space for cont_plot primitive.");
  whichobj->prim->which = CONT_PLOT;
  
  whichobj->prim->cont_plot = (contour_plot_type *)D_CALLOC(1,sizeof(contour_plot_type));
  if( !whichobj->prim->cont_plot )
    fatal("Can't get space for cont_plot.");


  if( !filename ){
    display("Look in the xterm...");
#ifndef USE_READLINE    
    printf("Enter the file name containing the contour plot data: ");
    scanf("%s",file_name);
#else
    theinline= readline("Enter the file name containing the contour plot data: ");
    add_history(theinline);
    if( theinline ){
      sscanf(theinline,"%s",file_name);
      free(theinline);
    } else {
      error("Bad file name");
      file_name[0] = 0;
    }
#endif
  } else{
    strcpy(file_name,filename);
  }
  

  /* open the file */
  infile = fopen(file_name,"r");
  if(!infile){
    printf("Problems opening file: %s\n",file_name);
    display("oooooops!");
    failed = 1;
  }

  if( !failed ){
    read_3D_data(infile,whichobj->prim->cont_plot);
    strcpy(whichobj->prim->cont_plot->filename,file_name);
  }

  /* check to see if any curves were actually read in.... */
  if(!whichobj->prim->cont_plot->num_curves || failed ){
    /* no... free the memory that we asked for */
    D_FREE(whichobj->prim->cont_plot);
    D_FREE(whichobj->prim);
    D_FREE(whichobj);
    whichobj=0;
    head->obj = 0;
    head = head->next;
  }
  else{

    /* set some defaults */
    whichobj->scale.x=whichobj->scale.y=whichobj->scale.z=2.0*GRAPHICS_SCALE;
    whichobj->cent.x=-180;whichobj->cent.y=180;
    whichobj->cent.z=0*GRAPHICS_SCALE;
    whichobj->trans.x=0;whichobj->trans.y=0;
    whichobj->trans.z=0;
    whichobj->prim->cont_plot->total_DOS_on_y = 1;
    whichobj->prim->cont_plot->num_levels = DEFAULT_NUM_OF_ZLEVELS;
    whichobj->prim->cont_plot->num_approx_pts = DEFAULT_NUM_APPROX_PTS;
    whichobj->prim->cont_plot->interp_kind = INTERP_CUBIC;
    whichobj->prim->cont_plot->order = DEFAULT_BSPLINE_ORDER;
    whichobj->prim->cont_plot->levels_kind = LEVELS_AUTO;
      
    /* make the button window */
#ifdef INTERACTIVE_USE    
    build_cont_plot_button_window(&button_wins,whichobj->prim->cont_plot);
    whichobj->prim->but_win = button_wins;
#endif

  }
  
}


/****************************************************************************
 *
 *                   Procedure draw_cont_plot
 *
 * Arguments: prim: pointer to primitive_type
 *             obj: pointer to object_type
 *            
 * Returns: none
 *
 * Action: Draws in a contout plot.
 *   
 *   if it's an FCO plot, the total and fragment DOSs are also plotted.
 *
 *   This takes care of the user specified scaling and translating.
 *
 ****************************************************************************/
void draw_cont_plot(prim_type *prim,object_type *obj)
{
  static XPoint *points=0;
  static point_type2D *fpoints=0;
  static int num_points=0;
  int i,j,num_2D_points;
  contour_plot_type *the_cont_plot;
  gnuplot_contour_type *contour;
  point_type2D origin,dim;
  char numstring[20];
  float xloc,yloc;
  float inv_xscale,inv_yscale;
  float xscale,yscale;
  int pts_so_far;
  float xref,yref;
  float xval,yval;
  float *max_x,*max_y,*min_x,*min_y,*max_z,*min_z;
  float *old_max_x,*old_max_y,*old_min_x,*old_min_y,*old_max_z,*old_min_z;

  if( !prim->cont_plot ) FATAL_BUG("Bogus primitive passed to draw_cont_plot.");
  the_cont_plot = prim->cont_plot;

  max_x = &(the_cont_plot->max_x);
  max_y = &(the_cont_plot->max_y);
  max_z = &(the_cont_plot->max_z);
  min_x = &(the_cont_plot->min_x);
  min_y = &(the_cont_plot->min_y);
  min_z = &(the_cont_plot->min_z);
  old_max_x = &(the_cont_plot->old_max_x);
  old_max_y = &(the_cont_plot->old_max_y);
  old_max_z = &(the_cont_plot->old_max_z);
  old_min_x = &(the_cont_plot->old_min_x);
  old_min_y = &(the_cont_plot->old_min_y);
  old_min_z = &(the_cont_plot->old_min_z);
  
  /* check to see if we need to re-determine the location of tic marks */
  if( *old_max_x != *max_x ||
     *old_max_y != *max_y ||
     *old_min_x != *min_x ||
     *old_min_y != *min_y ){


    preprocess_cont_plot_data(the_cont_plot);
    
    *old_max_x = *max_x;
    *old_max_y = *max_y;
    *old_max_z = *max_z;
    *old_min_x = *min_x;
    *old_min_y = *min_y;
    *old_min_z = *min_z;

  }
  
  /* check to see if we need memory for the Xpoints */
  if( the_cont_plot->num_x > the_cont_plot->num_y )
    num_2D_points = the_cont_plot->num_x;
  else
    num_2D_points = the_cont_plot->num_y;

  if( !points || num_points < num_2D_points ){
    if( points ) free(points);
    num_points = num_2D_points;
    points = (XPoint *)calloc(num_points+2,sizeof(XPoint));
    if( !points ) fatal("Can't allocate memory for point storage in draw_cont_plot.");
    if( fpoints ) free(fpoints);
    fpoints = (point_type2D *)calloc(num_points+2,sizeof(point_type2D));
    if( !fpoints ) fatal("Can't allocate memory for fpoint storage in draw_cont_plot.");    
  }
  
  /* determine the location of the origin on screen */
  origin.x = obj->cent.x + obj->trans.x + g_xmax / 2;
  origin.y = obj->cent.y - obj->trans.y + g_ymax / 2;  
  dim.x = DEF_CONT_PLOT_X * obj->scale.x;
  dim.y = DEF_CONT_PLOT_Y * obj->scale.y;
  
  /* scaling terms */
  inv_xscale = (the_cont_plot->max_x - the_cont_plot->min_x) / DEF_CONT_PLOT_X;
  inv_yscale = (the_cont_plot->max_y - the_cont_plot->min_y) / DEF_CONT_PLOT_Y;
  xscale = 1/inv_xscale;
  yscale = 1/inv_yscale;
  
  /* find the point in data space which will appear at the origin */
  xref = the_cont_plot->min_x * xscale;
  yref = the_cont_plot->min_y * yscale;

  /* determine the bounding box for this object */
  localmin.x = obj->bmin.x = origin.x;
  localmin.y = obj->bmin.y = origin.y - dim.y;
  localmax.x = obj->bmax.x = origin.x + dim.x;
  localmax.y = obj->bmax.y = origin.y;
  
  /* put in the tic marks now... */
  if( the_cont_plot->do_x_tics ){
    for(i=0;i<(int)rint(the_cont_plot->num_tics_x);i++){
      xloc = origin.x + obj->scale.x * (the_cont_plot->tic_start_x +
					i * the_cont_plot->tic_sep_x - xref);
      g_line(xloc,origin.y+obj->scale.y*TIC_DIM,xloc,origin.y);

      localmin.y += obj->scale.y*TIC_DIM;
      
      /*****
	
	do the labels
	
      *****/
      xval = (the_cont_plot->tic_start_x + i*the_cont_plot->tic_sep_x)*inv_xscale;
      if( fabs(xval) < 1e-12 ) xval = 0.0;
      sprintf(numstring,"%lg",xval);
      yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);

      g_center_text(xloc,yloc,numstring);
    }
  }

  /* Y tics */
  if( the_cont_plot->do_y_tics ){
    for(i=0;i<(int)rint(the_cont_plot->num_tics_y);i++){
      yloc = origin.y + obj->scale.y * (yref - the_cont_plot->tic_start_y -
					i * the_cont_plot->tic_sep_y);
      g_line(origin.x-obj->scale.x*TIC_DIM,yloc,origin.x,yloc);
      
      yval = (the_cont_plot->tic_start_y + i*the_cont_plot->tic_sep_y)*inv_yscale;
      if( fabs(yval) < 1e-12 ) yval = 0.0;
      sprintf(numstring,"%lg",yval);
      xloc = origin.x-obj->scale.x*TIC_DIM*1.1;

      g_right_text(xloc,yloc,numstring);
    }
  }  

  /* put in legends if they are needed */
  if( the_cont_plot->xlegend[0] != 0 && the_cont_plot->do_x_tics ){
    xloc = origin.x + dim.x/2;
    yloc = origin.y+obj->scale.y*TIC_DIM*(1.1);

    g_xlegend(xloc,yloc,the_cont_plot->xlegend);
  }

  if( the_cont_plot->ylegend[0] != 0 && the_cont_plot->do_y_tics ){
    yloc = origin.y - dim.y/2;
    xloc = origin.x-obj->scale.x*TIC_DIM*5.0;
    
    g_ylegend(xloc,yloc,the_cont_plot->ylegend);
  }

#if 1
  /* Now do the title */
  if( the_cont_plot->title[0] != 0 && the_cont_plot->do_title){
    xloc = origin.x + dim.x/2;
    yloc = origin.y - dim.y;

    g_title(xloc,yloc,the_cont_plot->title);
  }
#endif

  /* now do the plot itself */
  g_change_linewidth(0.25);
  for(i=0;i<the_cont_plot->num_curves;i++){

    if( the_cont_plot->curves_to_display[i] ){
      /* adjust the line styles as appropriate */
      if( the_cont_plot->type == CONT_FCO ){
	g_change_linestyle(the_cont_plot->styles[i]);
      }
      /* the contours are in a linked list, loop over that now */
      contour = the_cont_plot->contours[i];
      while(contour){
	if( the_cont_plot->type != CONT_FCO ){
	  if( contour->coords[0].z > 0 ){
	    g_change_linestyle(0);
	  } else if( contour->coords[0].z < 0 ){
	    g_change_linestyle(2);
	  } else{
	    g_change_linestyle(1);
	  }
	}
	  
	if( contour->coords[0].z <= the_cont_plot->max_z &&
	   contour->coords[0].z >= the_cont_plot->min_z ){
	   
	  /**********
	    
	    add points to the points array until we are either done with the contour 
	    or we've filled the array... then draw it.
	    
	    **********/
	  pts_so_far=0;
	  for(j=0;j<contour->num_pts;j++){
	    points[pts_so_far].x = fpoints[pts_so_far].x = origin.x +
	      (contour->coords[j].x - xref)*obj->scale.x;
	    points[pts_so_far].y = fpoints[pts_so_far].y = origin.y +
	      (yref - contour->coords[j].y)* obj->scale.y;
	    pts_so_far++;

	    if( pts_so_far == num_points ){
	      /* the array is full, draw it */
	      g_lines(points,fpoints,pts_so_far,0);
	      pts_so_far = 0;
	      /* put the current point back in the array */
	      points[pts_so_far].x = fpoints[pts_so_far].x = origin.x +
		(contour->coords[j].x - xref)*obj->scale.x;
	      points[pts_so_far].y = fpoints[pts_so_far].y = origin.y +
		(yref - contour->coords[j].y)* obj->scale.y;
	      pts_so_far++;
	    }
	  }
	  /* draw the points that are left over */
	  g_lines(points,fpoints,pts_so_far,0);
	  
	}
	/* move to the next contour */
	contour = contour->next;

      }
      /* set the line style back to the default value */
      g_change_linestyle(0);
    }
  }

  /* draw in a box surrounding the plot */
  g_change_linewidth(1);
  g_rectangle(origin.x,origin.y,dim.x,dim.y);

  /*******

    draw a line along the diagonal

  ********/
  if( the_cont_plot->type == CONT_FCO ){
    g_change_linestyle(2);
    if(*min_x == *min_y && *max_x == *max_y ){
      g_line(origin.x,origin.y,origin.x+dim.x,origin.y-dim.y);
    }
  }


  if( the_cont_plot->type == CONT_FCO ){
    /************
      
      okay, move over to the X edge of the plot and draw in the Total
      DOS
      
    *************/
    origin.x += dim.x;
    g_change_linestyle(0);
    pts_so_far = 0;

    for(j=0;j<the_cont_plot->raw_num_y;j++){
      xval = origin.x +
	the_cont_plot->data2D[j*(the_cont_plot->num_curves+1)].y*
	  obj->scale.x;
      yval = origin.y +
	(yref - the_cont_plot->data2D[j*(the_cont_plot->num_curves+1)].x)*
	  obj->scale.y;
      
      if( yval <= origin.y && yval >= origin.y - dim.y ){
	points[pts_so_far].x = fpoints[pts_so_far].x = xval;
	points[pts_so_far].y = fpoints[pts_so_far].y = yval;
	pts_so_far++;
	
	if( pts_so_far == num_points ){
	  /* the array is full, draw it */
	  g_lines(points,fpoints,pts_so_far,0);
	  pts_so_far = 0;
	  /* put the current point back in the array */
	  points[pts_so_far].x = fpoints[pts_so_far].x = xval;
	  points[pts_so_far].y = fpoints[pts_so_far].y = yval;
	  
	  pts_so_far++;
	}
      }
    }
    /* draw the points that are left over */
    g_lines(points,fpoints,pts_so_far,0);
    
    /************
      
      Do the same thing with the fragment DOS on the top edge
      
    *************/
    origin.x -= dim.x;
    origin.y -= dim.y;
    g_change_linestyle(0);
    pts_so_far = 0;
    for(i=1; i<=the_cont_plot->num_curves; i++){
      if( the_cont_plot->curves_to_display[i-1] ){
	g_change_linestyle(the_cont_plot->styles[i-1]);
	

	for(j=0;j<the_cont_plot->raw_num_x;j++){

	  xval = origin.x +
	    (the_cont_plot->data2D[j*(the_cont_plot->num_curves+1)+i].x-xref)*
	      obj->scale.x;
	  yval = origin.y +
	    (-the_cont_plot->data2D[j*(the_cont_plot->num_curves+1)+i].y)*
	      obj->scale.y;
	  if( xval >= origin.x && xval <= origin.x + dim.x ){
	    points[pts_so_far].x = fpoints[pts_so_far].x = xval;
	    points[pts_so_far].y = fpoints[pts_so_far].y = yval;
	    pts_so_far++;
      
	    if( pts_so_far == num_points ){
	      /* the array is full, draw it */
	      g_lines(points,fpoints,pts_so_far,0);
	      pts_so_far = 0;
	      /* put the current point back in the array */
	      points[pts_so_far].x = fpoints[pts_so_far].x = xval;
	      points[pts_so_far].y = fpoints[pts_so_far].y = yval;
	      pts_so_far++;
	    }
	  }
	}	
	/* draw the points that are left over */
	g_lines(points,fpoints,pts_so_far,0);
      }
    }
  }
  /* set the line style back to the default value */
  g_change_linestyle(0);
  

}

