#include <stdio.h>
#include <math.h>
#include <vector>
#include "shock_data_types.hpp"       /*tracer and shock data structers*/
#include "write_shock_catalogues.hpp" //correct catalogues
#include "timer.h"


tracer tin;                       /* buffer for adding tracers to tracer vectors */

struct interaction
{
  int snap_A;
  int snap_B;
  long idx_A;
  long id_A;
  long l_A;
  long o_A;
  float d_A;
  long idx_B;
  long id_B;
  long l_B;
  long o_B;
  float d_B;
  float frac_A;
  float frac_B;
  float frac_A_dense;
  float frac_B_dense;
  float x_A[3];
  float x_B[3];
  long n;
};
bool interaction_sort(interaction ia, interaction ib)
{
  //sort by decreasing frac_A
  return ia.frac_A>ib.frac_A;
}
void show_shocks(vector<shock> s);
void show_interactions(vector<interaction> ia);
void write_interactions(char fname[], vector<interaction> ia);
void keep_duplicates(vector<long> iunion, vector<long> *ioverlap);
void write_interaction_counts(int iA, int iB, long n_total, vector<long> n_ints, vector<long> o_ints);
int main(int argc, char **argv)
{
  char fdir[200];
  char fbase[200];
  char flist_A[200];
  char flist_B[200];
  char fdata_A[200];
  char fdata_B[200];
  char foutput[200];
  char fnints[200];

  int  imin = 795;
  int  imax = 801;

  int  iA;
  int  iB;

  long nA;
  long nB;

  long tt;
  long ss;

  float dr = sqrt(3)/512.;

  int nsearch = 10;

  float frac_A;
  float frac_B;


  float frac_A_dense;
  float frac_B_dense;
  float d_A;
  float d_B;
  float x_A[3];
  float x_B[3];

  vector<shock>  sA;
  vector<shock>  sB;
  vector<tracer> tA;
  vector<tracer> tB;

  vector<long> n_ints;
  vector<long> o_ints;
  long n_ints_total;
  long n_ints_snap;
  long n_ints_in;
  long n_ints_run=0;

  //vector<long> l_ia;
  //vector<long> o_ia;
  vector<interaction> ia;
  vector<interaction> ia_tmp;
  interaction ia_new;

  kdtree2 *bp_tree;
  kdtree2_result_vector res;
  array2dfloat bp_tree_data;
  vector<float> xc(3);

  float fdr = 3.; //how many cell distances?

  double time_A;
  double time_B;

  int *snap_list; //list of snapshots to use in tracking
  int n_snaps;    //total number of snapshots to use

  time_A = timer();

  //if we've supplied a range of snapshots
  //to search, use that
  if(argc>=3)
  {
  	imin = atoi(argv[1]);
  	imax = atoi(argv[2]);
  }

  if(argc>3)
  {
  	char fname_snap_list[200];



  	sprintf(fname_snap_list,"%s",argv[3]);

   	FILE *fp_snap_list;

   	if(!(fp_snap_list = fopen(fname_snap_list,"r")))
   	{
   		printf("Error opening %s\n",fname_snap_list);
   		exit(-1);
   	}

   	fscanf(fp_snap_list,"%d\n",&n_snaps);

    snap_list = (int *)malloc(n_snaps*sizeof(int));

   	for(int i=0;i<n_snaps;i++)
   	{
   		fscanf(fp_snap_list,"%d\n",&snap_list[i]);
   		//printf("snap_list[%d] = %d\n",i,snap_list[i]);
   	}
    fclose(fp_snap_list);

    imin = snap_list[n_snaps-1];
    imax = snap_list[0];
  }else{

    n_snaps = imax-imin+1;
    snap_list = (int *)malloc(n_snaps*sizeof(int));

    for(int i=0;i<n_snaps;i++)
      snap_list[i] = imax - i;

  }
  printf("imin = %d %d\n",imin,snap_list[n_snaps-1]);
  printf("imax = %d %d\n",imax,snap_list[0]);

  sprintf(fdir,"data/");
  sprintf(fbase,"peak.blended");
  sprintf(foutput,"interactions/interactions.%04d.%04d.txt",imin,imax);

  for(int isnap = 0; isnap<n_snaps-1; isnap++)
  {
  	iA = snap_list[isnap];
    iB = snap_list[isnap+1];

    printf("iA %04d\tiB %04d\n",iA,iB);

    sprintf(flist_A,"%s%s.%04d.list",fdir,fbase,iA);
    sprintf(flist_B,"%s%s.%04d.list",fdir,fbase,iB);	
    sprintf(fdata_A,"%s%s.%04d.dat", fdir,fbase,iA);
    sprintf(fdata_B,"%s%s.%04d.dat", fdir,fbase,iB);
    printf("%s\n%s\n%s\n%s\n",flist_A,flist_B,fdata_A,fdata_B);
    

    //if iA==imax, then read in the first
    //shock list from file.  If not, then
    //inherit the shock list from the last
    //shock list
  	if(iA==imax)
  	{
      read_shock_list(flist_A, &sA);
      nA = 0;
      for(size_t i=0;i<sA.size();i++)
		    nA += sA.size();
      tA.resize(nA);
      read_shock_data(fdata_A, sA, &tA);
      printf("sA->size() %ld\n",sA.size());
      printf("tA->size() %ld\n",tA.size());

  	}else{
  	  //copy sB into sA
  	  copy_shock_list(sB, &sA);
  	  copy_shock_data(tB, &tA);
  	  //free sB memory
  	  vector<shock>().swap(sB);
  	  vector<tracer>().swap(tB);  	  
  	}

  	printf("Reading shock data...\n");
 
  	//read sB list
    read_shock_list(flist_B, &sB);

    printf("sB->size() %ld\n",sB.size());

    //read sB data
    nB = 0;
    for(size_t i=0;i<sB.size();i++)
      nB += sB.size();
    tB.resize(nB);
    printf("tB.size() %ld\n",tB.size());
    read_shock_data(fdata_B, sB, &tB);

    printf("Done reading shocks...\n");

    //print out info about the shock list
    //show_shocks(sA);
    //show_shocks(sB);

    printf("tA %e\t%e\t%e\t%e\t%ld\n",tA[sA[0].o].x[0],tA[sA[0].o].x[1],tA[sA[0].o].x[2],tA[sA[0].o].d,tA[sA[0].o].id);
    printf("tB %e\t%e\t%e\t%e\t%ld\n",tB[sB[0].o].x[0],tB[sB[0].o].x[1],tB[sB[0].o].x[2],tB[sB[0].o].d,tB[sB[0].o].id);

    //fastest is probably to create a 
    //tree of the densest particles, and
    //then find the nearest peaks in B to A

    //construct tree for B peaks
    //we could construct this tree to account for
    //edge effects -- would require adding x6
    //or sorting to find edge cases
    bp_tree_data.resize(extents[sB.size()][3]);
    for(ss=0;ss<sB.size();ss++)
      for(int k=0;k<3;k++)
        bp_tree_data[ss][k] = tB[sB[ss].o].x[k];

    //build tree
    bp_tree = new kdtree2(bp_tree_data, true);

    //loop over A shocks
    for(long ss=0;ss<sA.size();ss++)
    //for(ss=0;ss<100;ss++)
    {
      n_ints_in = 0;
      if(sA[ss].l>1)
      {
        //this shock might have interactions
        long nkeep=100;

        //compute the density location
        d_A = 0;
        for(int k=0;k<3;k++)
          x_A[k] = 0;
        for(long tt=0;tt<sA[ss].l;tt++)
        {
          if(tA[sA[ss].o+tt].d>0.9*tA[sA[ss].o].d && tt>=nkeep)
            break;
          d_A += tA[sA[ss].o+tt].d;
          for(int k=0;k<3;k++)
            x_A[k] += tA[sA[ss].o+tt].d*tA[sA[ss].o+tt].x[k];
        }
        for(int k=0;k<3;k++)
          x_A[k] /= d_A;

        //printf("******\n");

        //perform a search on the location
        //of the A peak
        for(int k=0;k<3;k++)
          xc[k] = tA[sA[ss].o].x[k];

	    //find the 10 closest peaks
	    //how to deal with wrapping?
	    //we could check proximity to edges
	    //and then wrap if necessary

	    //do search
	    //bp_tree->n_nearest(xc,10,res);
	    bp_tree->r_nearest(xc,fdr*dr*dr,res);
      int force_flag = 0;
      int n_force = 5;
      if(res.size()<n_force)
      {
        bp_tree->n_nearest(xc,n_force,res);
        force_flag = 1;
      }

	    //print the results
	    //printf("***** iA %ld l %ld res.size() %ld\n",ss,sA[ss].l,res.size());
      //if(force_flag)
        //printf("FORCED\n");
	    //for(size_t i=0;i<res.size();i++)
	  	  //printf("A %10ld B %10ld Ad %5.4e Bd %5.4e dis %5.4e dr %5.4e\n",sA[ss].id,tB[sB[res[i].idx].o].id,sA[ss].d,tB[sB[res[i].idx].o].d,sqrt(res[i].dis),dr);

	    //keep track of the overlap
	    vector<long> idA;
      vector<long> idAdense;
      long n_dense_min = 100; 
	    for(tt=0;tt<sA[ss].l;tt++)
      {
	  	  idA.push_back(tA[sA[ss].o+tt].id);
        if(tA[sA[ss].o+tt].d > tA[sA[ss].o].d*0.9 || tt<n_dense_min)
          idAdense.push_back(tA[sA[ss].o+tt].id);
      }

   	    //sort idA
	    std::sort(idA.begin(),idA.end());
      std::sort(idAdense.begin(),idAdense.end());

	    //loop over results
	    for(size_t i=0;i<res.size();i++)
	    {
	  	  //if peak moved by less than dr
	  	  if(1)
	  	  //if(res[i].dis<fdr*dr*dr)
	  	  {
	  	    //compute the density location
	  	    d_B = 0;
	  	    for(int k=0;k<3;k++)
	  	      x_B[k] = 0;
	  	    for(long tt=0;tt<sB[res[i].idx].l;tt++)
	  	    {
	  		  if(tB[sB[res[i].idx].o+tt].d>0.9*tB[sB[res[i].idx].o].d && tt>=nkeep)
	  			  break;
	  		  d_B += tB[sB[res[i].idx].o+tt].d;
	  		  for(int k=0;k<3;k++)
	  		    x_B[k] += tB[sB[res[i].idx].o+tt].d*tB[sB[res[i].idx].o+tt].x[k];
	  	    }
	  	    for(int k=0;k<3;k++)
	  	      x_B[k] /= d_B;

	  	    //record the ids from sB
	  	    vector<long> idB, io, idBdense, iodense;
		      for(tt=0;tt<sB[res[i].idx].l;tt++)
          {
		        idB.push_back(tB[sB[res[i].idx].o+tt].id);
            if(tB[sB[res[i].idx].o+tt].d > tA[sB[res[i].idx].o].d*0.9 || tt<n_dense_min)
              idBdense.push_back(tB[sB[res[i].idx].o+tt].id);
          }
		      //add the ids from A
		      for(tt=0;tt<idA.size();tt++)
		  	   idB.push_back(idA[tt]);
          for(tt=0;tt<idAdense.size();tt++)
           idBdense.push_back(idAdense[tt]);

		      //sort union
		      std::sort(idB.begin(),idB.end());
          std::sort(idBdense.begin(),idBdense.end());


		      //keep only the duplicates
		      keep_duplicates(idB,&io);
          keep_duplicates(idBdense,&iodense);

		    //compare the ids, check to see how many are duplicated
		    frac_A = ((double) io.size())/((double) sA[ss].l);
		    frac_B = ((double) io.size())/((double) sB[res[i].idx].l);
        frac_A_dense = ((double) iodense.size())/((double) idAdense.size());
        frac_B_dense = ((double) iodense.size())/((double) idBdense.size());
		    if(io.size()>0)
		    {
	  	      //this is a potential peak
	  	      //printf("A %10ld B %10ld Ad %5.4e Bd %5.4e dis %5.4e dr %5.4e\n",sA[ss].id,tB[sB[res[i].idx].o].id,sA[ss].d,tB[sB[res[i].idx].o].d,sqrt(res[i].dis),sqrt(fdr)*dr);

           //printf("iosize %ld idBsize %ld idAsize %ld iodsize %ld idBdsize %ld idAdsize %ld\n",io.size(),idB.size(),idA.size(),iodense.size(),idBdense.size(),idAdense.size());
  		     //printf("ia %10ld\tib %10ld\tida %10ld\t idb %10ld\tio %10ld\tfrac A %5.4e\tfrac B %5.4e frac Ad %5.4e\tfrac Bd %5.4e\n",sA[ss].l,sB[res[i].idx].l,sA[ss].id,sB[res[i].idx].id,io.size(),frac_A,frac_B,frac_A_dense,frac_B_dense);
  		     //record the interaction
  		     ia_new.snap_A = iA;
  		     ia_new.snap_B = iB;
  		     ia_new.idx_A  = ss;
  		     ia_new.id_A   = sA[ss].id;
  		     ia_new.d_A    = sA[ss].d;
  		     ia_new.l_A    = sA[ss].l;
  		     ia_new.o_A    = sA[ss].o;
  		     ia_new.idx_B  = res[i].idx;
  		     ia_new.id_B   = sB[res[i].idx].id;
  		     ia_new.d_B    = sB[res[i].idx].d;
  		     ia_new.l_B    = sB[res[i].idx].l;
  		     ia_new.o_B    = sB[res[i].idx].o;
  		     ia_new.n      = io.size();
  		     ia_new.frac_A = frac_A;
  		     ia_new.frac_B = frac_B;
           ia_new.frac_A_dense = frac_A_dense;
           ia_new.frac_B_dense = frac_B_dense;
  		     for(int k=0;k<3;k++)
  		     {
  		     	ia_new.x_A[k] = x_A[k];
  		     	ia_new.x_B[k] = x_B[k];
  		     }
  		     ia_tmp.push_back(ia_new);
		    }
		    //sort the temporary interactions by decreasing frac_A
		    std::sort(ia_tmp.begin(),ia_tmp.end(),interaction_sort);

		    //store these interactions
		    for(int i=0;i<ia_tmp.size();i++)
		  	  ia.push_back(ia_tmp[i]);


		  	n_ints_in += ia_tmp.size();

		    //free the temporary interactions
		    vector<interaction>().swap(ia_tmp);
	  	    //free the ids
	  	    vector<long>().swap(idB);
          vector<long>().swap(idBdense);  	  
	  	  }//end if of valid interaction
	    }//end loop over search results

	    //free the ids
	    vector<long>().swap(idA);
      vector<long>().swap(idAdense);
	   }else{ //end if s.l>1
      n_ints_in = 0; //explicitly not interested in l=1 peaks
     }

	   //keep track of how many
	   //interactions this shock had
	   n_ints.push_back(n_ints_in);
     n_ints_run += n_ints_in;
     //printf("iA %ld n_ints_in %ld n_ints_run %ld sA[ss].l %ld\n",ss,n_ints_in,n_ints_run,sA[ss].l);

     //if(n_ints_in==0 && sA[ss].l>1)
       //exit(-1);
    }//end loop over sA

    //free tree and associated data
    bp_tree_data.resize(extents[0][0]);
    free(bp_tree);

    printf("sA.size() %ld n_ints.size() %ld\n",sA.size(),n_ints.size());

  	//keep track of the offsets for each
  	//shock in the list of interactions
    o_ints.resize(n_ints.size());
    o_ints[0] = 0;
    n_ints_total = n_ints[0];
    for(int i=1;i<n_ints.size();i++)
    {
      n_ints_total += n_ints[i];
  	  o_ints[i]     = o_ints[i-1] + n_ints[i-1];
  	}
    printf("n_ints_total %ld, n_ints_run %ld, n_ints.size() %ld, sA.size() %ld\n",n_ints_total,n_ints_run,n_ints.size(),sA.size());
  	write_interaction_counts(iB,iA,n_ints_total,n_ints,o_ints);

	vector<long>().swap(n_ints);
	vector<long>().swap(o_ints);
  }//end loop over snapshots



  //for(int i=0;i<ia.size();i++)
	//printf("snap_A %04d\tia %10ld\tib %10ld\tn %10ld\tfrac A %5.4e\tfrac B %5.4e\txA %e %e %e\txB %e %e %e\n",ia[i].snap_A,ia[i].id_A,ia[i].id_B,ia[i].n,ia[i].frac_A,ia[i].frac_B,ia[i].x_A[0],ia[i].x_A[1],ia[i].x_A[2],ia[i].x_B[0],ia[i].x_B[1],ia[i].x_B[2]);

  //save the interactions to a file
  write_interactions(foutput,ia);

  time_B = timer();

  printf("Total time = %es.\n",time_B-time_A);

  return 0;
}
void show_shocks(vector<shock> s)
{
  int nlim = 10;
  printf("***************\n");
  for(int i=0;i<nlim;i++)
	  printf("i %6d\tl %10ld\to %10ld\td %15.14e\tid %10ld\tmin %5.4e %5.4e %5.4e\tmax %5.4e %5.4e %5.4e\n",i,s[i].l,s[i].o,s[i].d,s[i].id,s[i].min[0],s[i].min[1],s[i].min[2],s[i].max[0],s[i].max[1],s[i].max[2]);
  printf("***************\n");
  for(int i=s.size()-nlim;i<s.size();i++)
	  printf("i %6d\tl %10ld\to %10ld\td %15.14e\tid %10ld\tmin %5.4e %5.4e %5.4e\tmax %5.4e %5.4e %5.4e\n",i,s[i].l,s[i].o,s[i].d,s[i].id,s[i].min[0],s[i].min[1],s[i].min[2],s[i].max[0],s[i].max[1],s[i].max[2]);

}
void show_interactions(vector<interaction> ia, vector<long> n_ints, vector<long> o_ints)
{
	int nlim = ia.size();
	printf("**********\n");
	for(int i=0;i<nlim;i++)
	{
		printf("snap_A %04d\tsnap_B %04d\tidA %10ld\tidB %10lddA %5.4e\tdB %5.4e\n",ia[i].snap_A,ia[i].snap_B,ia[i].id_A,ia[i].id_B,ia[i].d_A,ia[i].d_B);
	}
}
void keep_duplicates(vector<long> iunion, vector<long> *ioverlap)
{
  vector<long>::iterator ia;

  ia = std::adjacent_find(iunion.begin(), iunion.end());
  if(ia!=iunion.end())
  {
    ioverlap->push_back(*ia);
    while(ia!=iunion.end())
    {
      ia = std::adjacent_find(++ia, iunion.end());
      if(ia!=iunion.end())
        ioverlap->push_back(*ia);
    }
  }
}
void write_interactions(char fname[], vector<interaction> ia)
{
  FILE *fp;
  if(!(fp=fopen(fname,"w")))
  {
  	printf("Error opening %s.\n",fname);
  	exit(-1);
  }
  fprintf(fp,"%ld\n",ia.size());
  for(size_t i=0;i<ia.size();i++)
  {
  	fprintf(fp,"%04d %04d %8ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e %10ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e %10ld %8ld %8ld %5.4e %5.4e %5.4e %5.4e\n",ia[i].snap_A,ia[i].snap_B,ia[i].idx_A,ia[i].idx_B,ia[i].n,ia[i].frac_A,ia[i].frac_B,ia[i].frac_A_dense,ia[i].frac_B_dense,ia[i].id_A,ia[i].l_A,ia[i].o_A,ia[i].d_A,ia[i].x_A[0],ia[i].x_A[1],ia[i].x_A[2],ia[i].id_B,ia[i].l_B,ia[i].o_B,ia[i].d_B,ia[i].x_B[0],ia[i].x_B[1],ia[i].x_B[2]);
  }
  fclose(fp);
}
void write_interaction_counts(int iA, int iB, long n_total, vector<long> n_ints, vector<long> o_ints)
{
  FILE *fp;
  char fname[200];

  sprintf(fname,"interactions/interaction_count.%04d.%04d.txt",iA,iB);
  printf("wic: n_ints.size() %ld\n",n_ints.size());
  if(!(fp=fopen(fname,"w")))
  {
  	printf("Error opening %s.\n",fname);
  	exit(-1);
  }
  //fprintf(fp,"%ld\n",n_total);
  fprintf(fp,"%ld\n",n_ints.size());
  for(int i=0;i<n_total;i++)
  	fprintf(fp,"%ld\t%ld\n",n_ints[i],o_ints[i]);
  fclose(fp);
}
