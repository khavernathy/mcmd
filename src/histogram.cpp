#include <sys/types.h>
#ifdef WINDOWS
	#include <io.h>
#else
	#include <unistd.h>
#endif


void output(char *msg) {
    printf("%s", msg);
    fflush(stdout);
}

/* Take a vector in cartesian coordinates (cart[]) and
 * convert it to fractional coordinates.
 * store the answer in anwer[]     */
int cart2frac(double *answer, double *cart, System &system){

        int i,j;

        for(i=0;i<3;i++) answer[i]=0.0;

        for(i=0;i<3;i++){
                for(j=0;j<3;j++){
		//we use transpose(recip_basis), because transpose(recip_basis).basis_vector = <1,0,0> , <0,1,0> or <0,0,1>
                        answer[i]+=system.pbc.reciprocal_basis[j][i]*cart[j];
                }
        }

        return 1;
}


/* Take a vector in fractional coordinates (frac[]) and
 * convert it to cartesian coordinates.
 * store the answer in answer[]     */
int frac2cart(double *answer, double *frac, System &system)
{

	int i,j;

	for(i=0;i<3;i++) answer[i]=0.0;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			answer[i]+=system.pbc.basis[j][i]*frac[j];
		}
	}


	return 1;
}



/* update the avg_histogram stored on root */
void update_root_histogram(System &system)
{
	int i,j,k;
	int xdim=system.grids.histogram->x_dim;
	int ydim=system.grids.histogram->y_dim;
	int zdim=system.grids.histogram->z_dim;

	for(k=0; k < zdim; k++){
		for(j=0; j < ydim; j++){
			for(i=0; i < xdim; i++){
				system.grids.avg_histogram->grid[i][j][k]+=system.grids.histogram->grid[i][j][k];
				/* norm_total is updated here to normalize upon write out */
				system.grids.avg_histogram->norm_total+=system.grids.histogram->grid[i][j][k];
			}
		}
	}
}


/* zero out the histogram grid */
void zero_grid(int ***grid, System &system)
{
	int i,j,k;
	int xdim=system.grids.histogram->x_dim;
	int ydim=system.grids.histogram->y_dim;
	int zdim=system.grids.histogram->z_dim;

	for(k=0; k < zdim; k++){
		for(j=0; j < ydim; j++){
			for(i=0; i < xdim; i++){
				grid[i][j][k]=0;
			}
		}
	}
}
/* compute the histogram bin number to place this molecule in */
void compute_bin(double *cart_coords, System &system, int *bin_vector)
{
	double frac_coords[3];
	int Abin,Bbin,Cbin;
	/* we need the fractional coords to simplify the bin calculation */
	cart2frac(frac_coords,cart_coords,system);

	/* the coordinate system in the simulation is from -0.5 to 0.5 */
	/* so we need to correct for this: add 0.5 to each dimension */
	frac_coords[0]+=0.5;
	frac_coords[1]+=0.5;
	frac_coords[2]+=0.5;

	/* compute bin in each dimension */
	Abin=(int)floor(frac_coords[0]*system.grids.histogram->x_dim);
	Bbin=(int)floor(frac_coords[1]*system.grids.histogram->y_dim);
	Cbin=(int)floor(frac_coords[2]*system.grids.histogram->z_dim);

	/* return result to the bin_vector passed in */
	bin_vector[0]=Abin;
	bin_vector[1]=Bbin;
	bin_vector[2]=Cbin;
}

void wrap1coord(double *unwrapped, double *wrapped, System &system)
{
	double unit[3],offset[3],frac[3];
	int i,j;

	/* zero the vectors */
	for(i=0;i<3;i++){ offset[i]=0; unit[i]=0; frac[i]=0;}

	/* put coords in fractional representation */
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			//we use transpose(recip_basis), because transpose(recip_basis).basis_vector = <1,0,0> , <0,1,0> or <0,0,1>
			frac[i]+=system.pbc.reciprocal_basis[j][i]*unwrapped[j];

	/* any fractional coord > .5 or < -.5 round to 1,-1 etc. */
	for(i=0;i<3;i++) unit[i] = rint(frac[i]);

	/* multiply this rounded fractional unit vector by basis */
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			offset[i]+=system.pbc.basis[j][i]*unit[j];

	/* subtract this distance from the incoming vector */
	for(i=0;i<3;i++) wrapped[i] = unwrapped[i] - offset[i];

}

/* population histogram should be performed only every corr time
 * only root will compile the sum into grids.avg_histogram
 * the node only grid.histogram will be rezeroed at every corr time
 * to prevent overflow.
 * NOTE: The histogram is normalized to 1.
 * e.g. if there are a total of 328 bin counts . every bin is
 * divided by 328. */
void population_histogram(System &system)
{
	int i;
	int bin[3];
	double wrappedcoords[3];
	for(i=0; i<system.molecules.size(); i++){
		if(!system.molecules[i].frozen){
			/* wrap the coordinates of mol_p */
			wrap1coord(system.molecules[i].com,wrappedcoords,system);
			/* compute what bin to increment. store answer in bin[] */
			compute_bin(wrappedcoords,system,bin);
			/* increment the bin returned in bin[] */
			(system.grids.histogram->grid[(bin[0])][(bin[1])][(bin[2])])++;
		}
	}

}
/* This writes out the grid with the Cbin (last index) varying
 * the fastest.  Line break between dimensions, double line
 * break between ZY sheets, ####'s between complete sets
 * Remember, for this to work, we had to offset the origin
 * by 1/2 a bin width */
void write_histogram(FILE *fp_out, int ***grid, System &system)
{
	int i,j,k;
	int xdim,ydim,zdim;
	int count=0;

	xdim=system.grids.histogram->x_dim;
	ydim=system.grids.histogram->y_dim;
	zdim=system.grids.histogram->z_dim;

	rewind(fp_out);
	fprintf(fp_out,"# OpenDX format population histogram\n");
	fprintf(fp_out,"object 1 class gridpositions counts %d %d %d\n",xdim,ydim,zdim);
	fprintf(fp_out,"origin\t%f\t%f\t%f\n",system.grids.histogram->origin[0],system.grids.histogram->origin[1],system.grids.histogram->origin[2]);
	fprintf(fp_out,"delta \t%f\t%f\t%f\n",system.grids.histogram->delta[0][0],system.grids.histogram->delta[0][1],system.grids.histogram->delta[0][2]);
	fprintf(fp_out,"delta \t%f\t%f\t%f\n",system.grids.histogram->delta[1][0],system.grids.histogram->delta[1][1],system.grids.histogram->delta[1][2]);
	fprintf(fp_out,"delta \t%f\t%f\t%f\n",system.grids.histogram->delta[2][0],system.grids.histogram->delta[2][1],system.grids.histogram->delta[2][2]);
	fprintf(fp_out,"\n");
	fprintf(fp_out,"object 2 class gridconnections counts %d %d %d\n",xdim,ydim,zdim);
	fprintf(fp_out,"\n");

	for(i=0; i < xdim; i++){
		for(j=0; j < ydim; j++){
			for(k=0; k < zdim; k++){
				fprintf(fp_out,"%f ",(float)(grid[i][j][k]) / (float)(system.grids.avg_histogram->norm_total));
				count+=grid[i][j][k];
			}
			fprintf(fp_out,"\n");
		}
		fprintf(fp_out,"\n");
	}

	fprintf(fp_out,"# count=%d\n",count);
	fprintf(fp_out,"attribute \"dep\" string \"positions\"\n");
	fprintf(fp_out,"object \"regular positions regular connections\" class field\n");
	fprintf(fp_out,"component \"positions\" value 1\n");
	fprintf(fp_out,"component \"connections\" value 2\n");
	fprintf(fp_out,"component \"data\" value 3\n");
	fprintf(fp_out,"\nend\n");
	fflush(fp_out);
}

/* Returns the magnitude of a 3-vector */
double magnitude(double *vector){
	int i;
	double mag=0;
	for(i=0;i<3;i++) mag+=vector[i]*vector[i];
	return sqrt(mag);
}

void allocate_histogram_grid(System &system)
{
	int i,j;
	int x_dim=system.grids.histogram->x_dim;
	int y_dim=system.grids.histogram->y_dim;
	int z_dim=system.grids.histogram->z_dim;

	/* allocate a 3D grid for the histogram */
	system.grids.histogram->grid = (int ***) calloc(x_dim, sizeof(int **));
	for(i=0; i<x_dim; i++) {
		system.grids.histogram->grid[i] = (int **) calloc(y_dim, sizeof(int *));
	}


	for(i=0; i<x_dim; i++){
		for(j=0; j<y_dim; j++){
			system.grids.histogram->grid[i][j] = (int *) calloc(z_dim, sizeof(int));
		}
	}

	/* allocate an avg_histogram grid */
		system.grids.avg_histogram->grid = (int ***) calloc(x_dim, sizeof(int **));
		for(i=0; i<x_dim; i++) {
			system.grids.avg_histogram->grid[i] = (int **) calloc(y_dim, sizeof(int *));
		}
		for(i=0; i<x_dim; i++){
			for(j=0; j<y_dim; j++){
				system.grids.avg_histogram->grid[i][j] = (int *) calloc(z_dim, sizeof(int));
			}
		}


}

/* These are needed by dxwrite routines.
 * delta is a transformation matrix that defines the
 * step size for setting up a grid in OpenDX
 * see the DX userguide.pdf appendix B for
 * more details. */
void setup_deltas(histogram_t *hist, System &system)
{
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			hist->delta[i][j]=system.pbc.basis[j][i]/hist->count[i];
	/* we divide by the count to get our actual step size in each dimension */

    /* dfranz: IMPORTANT: CORRECTIONS FOR NON-ORTHORHOMBIC SYSTEMS */
    // that is to say, the above is fine for cubic systems but results 
    // are a weird and rotated DX cell for non 90/90/90
    hist->delta[1][0] = hist->delta[0][1];
        hist->delta[0][1] = 0.0;
    hist->delta[2][0] = hist->delta[0][2];
        hist->delta[0][2] = 0.0;

}

/* Because OpenDX puts our data values at the points of the grid,
 * we need to offset them by half a bin width to reflect the fact that
 * we have a histogram. Our data really lies between each point. */
void offset_dx_origin(double *real_origin_cartesian, histogram_t *hist, System &system)
{
	double fractional_binwidth[3];
	double cart_halfbin[3];

	/* figure out how wide each bin is in each dimension */
	fractional_binwidth[0] = 1.0 / hist->x_dim;
	fractional_binwidth[1] = 1.0 / hist->y_dim;
	fractional_binwidth[2] = 1.0 / hist->z_dim;

	/* figure out how wide half a binwidth is in cartesians */
	fractional_binwidth[0]/=2.0;
	fractional_binwidth[1]/=2.0;
	fractional_binwidth[2]/=2.0;
	frac2cart(cart_halfbin, fractional_binwidth, system);

	/* add this value to the origin */
	real_origin_cartesian[0]+=cart_halfbin[0];
	real_origin_cartesian[1]+=cart_halfbin[1];
	real_origin_cartesian[2]+=cart_halfbin[2];

}

/* Variables needed upon printing in openDX native format.
 * see DX userguide.pdf appendix B for details. */
void setup_dx_variables(histogram_t *hist, System &system)
{
	double vec[3],origin[3];

	/* setup counts */
	hist->count[0]=hist->x_dim;
	hist->count[1]=hist->y_dim;
	hist->count[2]=hist->z_dim;

	/* setup origin */
	vec[0]=-0.5; vec[1]=-0.5; vec[2]=-0.5;
	frac2cart(origin, vec, system);

	/* IMPORTANT!!! I am offsetting the origin by 1/2 a bin in each dimension */
	/* the result of origin is not the true origin!!! */
	offset_dx_origin(origin,hist,system);

	hist->origin[0]=origin[0];
	hist->origin[1]=origin[1];
	hist->origin[2]=origin[2];

	/* setup deltas */
	setup_deltas(system.grids.histogram,system);

	/* setup N data points */
	hist->n_data_points = hist->x_dim * hist->y_dim * hist->z_dim;
}

/* Setup the various quantities that define the histogram grid */
void setup_histogram(System &system){

	char linebuf[256];
	double trial_vec1[3];
	double trial_vec2[3];
	double magA,magB,magC;
	int Nbins,x_dim,y_dim,z_dim;

	/* get the magnitudes of all the basis vectors and test the frac2cart routine.
	 * define a fractional vector (1,0,0) and transform it with our basis.
	 * then calculate its magnitude.  Do this in all 3 dimensions         */
	trial_vec1[0]=1.0;	trial_vec1[1]=0.0;	trial_vec1[2]=0.0;
	frac2cart(trial_vec2,trial_vec1,system);
	magA = magnitude(trial_vec2);

	trial_vec1[0]=0.0;	trial_vec1[1]=1.0;	trial_vec1[2]=0.0;
	frac2cart(trial_vec2,trial_vec1,system);
	magB = magnitude(trial_vec2);

	trial_vec1[0]=0.0;	trial_vec1[1]=0.0;	trial_vec1[2]=1.0;
	frac2cart(trial_vec2,trial_vec1,system);
	magC = magnitude(trial_vec2);

	/* calculate the number of bins in each fractional coordinate */
	x_dim=rint(magA/system.hist_resolution);
	y_dim=rint(magB/system.hist_resolution);
	z_dim=rint(magC/system.hist_resolution);
	sprintf(linebuf,"HISTOGRAM: %f A; resolution -> %d bins(A) * %d bins(B) * %d bins(C)\n",system.hist_resolution,x_dim,y_dim,z_dim);output(linebuf);

	Nbins=x_dim * y_dim * z_dim;
	sprintf(linebuf,"HISTOGRAM: Total Bins = %d\n",Nbins); output(linebuf);


	system.grids.histogram->x_dim = x_dim;
	system.grids.histogram->y_dim = y_dim;
	system.grids.histogram->z_dim = z_dim;
	system.grids.avg_histogram->x_dim = x_dim;
	system.grids.avg_histogram->y_dim = y_dim;
	system.grids.avg_histogram->z_dim = z_dim;
	system.grids.histogram->n_data_points = Nbins;
	system.n_histogram_bins = Nbins;
	system.grids.avg_histogram->norm_total=0;
	system.grids.histogram->norm_total=0;

	setup_dx_variables(system.grids.histogram,system);
}
