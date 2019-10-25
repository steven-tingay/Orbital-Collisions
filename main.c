#include "main.h"
int main (int argc, char *argv[])
{
	int suc,yesNo,loops;
	long double charLLen,mass,totalMass;
	char *stopString;
	srand(time(NULL));
	if(argc==1)
	{	
		printf("\nChoose a lower bound characteristic Length(m)\n\n");
		suc=scanf("%Lf",&charLLen);
		printf("\nChoose an impact mass(kg)\n\n");
		suc=scanf("%Lf",&mass);
		printf("\nWould you like all fragments to be outputed to terminal? 0=yes,1=no\n");
		suc=scanf("%d",&yesNo);
	}
	else
	{
		charLLen = strtold(argv[1],&stopString);
		mass = strtod(argv[2],&stopString);
		yesNo=atoi(argv[3]);
		loops = atoi(argv[4]);	
	}	
	if(loops<=1)
	{
		singleSim(charLLen, mass,yesNo);
	}
	else
	{
		multiSim(charLLen,mass, loops);
	}		
	return 0;
}
void singleSim(long double charLLen, long double mass, int yesNo)
{
	linkedList** list;
	long double totalMass;
	list = createLL(charLLen,mass);
	totalMass=initiateAM(list,1);
	assignDelV(list,1);
	if(yesNo==0)
	{
		print(list);
	}
	printFile(list);
	deAllocLL(list);
	printf("Total Mass = %Lf\n",totalMass);
}	
void multiSim(long double charLLen, long double mass,int loops)
{
	linkedList **list;
	long double totalMass=0,previous=0.0;
	int ii=0;
	while(ii<loops)
	{
		list = createLL(charLLen,mass);
		totalMass = initiateAM(list,1);
		totalMass-=previous;
		previous+= totalMass;
		assignDelV(list,1);
		printf("Total Mass = %Lf\n",totalMass);
		deAllocLL(list);
		ii++;
	}
}	
long double initiateAM(linkedList** list,int impactType)
{
	int ii=0,numBins;
	numBins=list[0]->numbBins;
	long double ** probDF,** cumulDF,parArray[5],totalMass;
	probDF = (long double**)calloc(2,sizeof(long double*));
	createDArray(2,MAX_PDF,probDF);
	setStanVal(probDF,MAX_PDF,X_LOWER,X_UPPER);
	cumulDF=(long double**)calloc(2,sizeof(long double*));
	createDArray(2,MAX_PDF/4,cumulDF);
	setStanVal(cumulDF,MAX_PDF/4,X_LOWER+((X_UPPER-X_LOWER)/(MAX_PDF/4)),X_UPPER);
	while(ii<numBins)
	{
		totalMass=assignAM(list,ii,impactType,probDF,cumulDF);
		ii++;
	}
	deAllocA(probDF,2);
	deAllocA(cumulDF,2);
	return totalMass;
}
long double assignAM(linkedList** list,int arrayIndex,int impactType,long double** probDF,long double** cumulDF)
{
	int ii,numNodes;
	static int counter=0;
	static long double totalMass=0;
	long double parArray[5];
	numNodes = list[arrayIndex]->numbNodes;
	if(impactType==0)
	{
		getRBParam(log10(list[arrayIndex]->binCLen),parArray);
	}
	else
	{
		getSCParam(log10(list[arrayIndex]->binCLen),parArray);
	}
	addVal(probDF,parArray,MAX_PDF);
	createCDF(probDF, cumulDF, MAX_PDF,X_LOWER,X_UPPER);
	for(ii=0;ii<numNodes;ii++)
	{
		list[arrayIndex]->head[ii]->areaMass=sampleCDF(cumulDF);
		list[arrayIndex]->head[ii]->charLen=list[arrayIndex]->binCLen;
		list[arrayIndex]->head[ii]->csArea = 0.556945*pow(list[arrayIndex]->binCLen,2.0047077);
		list[arrayIndex]->head[ii]->mass = list[arrayIndex]->head[ii]->csArea/
			list[arrayIndex]->head[ii]->areaMass;
		totalMass+=list[arrayIndex]->head[ii]->mass;
	}
	//printf("cumulative bin[%d] mass = %Lf\n",counter,totalMass);
	counter++;
	return totalMass;
}
void assignDelV(linkedList** list, int impactType)
{
	int ii=0,aa=0,numBins,numNodes;
	long double **probDF,**cumulDF;
	probDF = (long double**)calloc(2,sizeof(long double*));
	cumulDF= (long double**)calloc(2,sizeof(long double*));
	createDArray(2,MAX_PDF,probDF);
	createDArray(2,MAX_PDF/4,cumulDF);
	setStanVal(cumulDF,MAX_PDF/4,XV_LOWER+((XV_UPPER-XV_LOWER)/(MAX_PDF/4)),XV_UPPER);
	numBins=list[0]->numbBins;
	while(ii<numBins)
	{
		numNodes=list[ii]->numbNodes;
		while(aa<numNodes)
		{
			addVelVal(probDF,log10(list[ii]->head[aa]->areaMass),MAX_PDF,impactType);
			createCDF(probDF,cumulDF,MAX_PDF,XV_LOWER,XV_UPPER);
			list[ii]->head[aa]->deltaV=sampleCDF(cumulDF);
			aa++;
		}
		aa=0;
		ii++;
	}
	deAllocA(probDF,2);
	deAllocA(cumulDF,2);
}
linkedList** createLL(long double charLLen,double mass)
{
	int numNodes,numBins=0,ii=0,aa=0;
	long double newUpper;
	linkedList** list;
	fragment* newNode;
	newUpper = charLenUp(charLLen,mass);
	numBins=(int)((newUpper-charLLen)*100);
	printf("numBins = %d, Upper char length = %Lf m\n",numBins,newUpper);
	list=(linkedList**)malloc(numBins*sizeof(linkedList*));
	while(charLLen<newUpper)
	{
		list[ii]=(linkedList*)malloc(sizeof(linkedList));
		numNodes = round(0.1*pow(mass,0.75)*(pow(charLLen,-1.71)-pow(charLLen+0.01,-1.71)));
		list[ii]->head=(fragment**)malloc(numNodes*sizeof(fragment*));
		while(aa<numNodes)
		{
			list[ii]->head[aa]=(fragment*)malloc(sizeof(fragment));
			aa++;
		}
		aa=0;
		list[ii]->numbNodes= numNodes; 
		list[ii]->numbBins = numBins;
		list[ii]->binCLen = charLLen;
		charLLen+=0.01;
		ii++;
	}
	return list;
}
void print(linkedList** list)
{
	int ii=0,aa=0,numNodes,numBins;
	numBins=list[ii]->numbBins;
	while(ii<numBins)
	{
		numNodes=list[ii]->numbNodes;
		printf("#############\n\n Bin:%d\n\n############\n",ii); 
		while(aa<numNodes)
		{
			printf("**************************************\n\n Fragment %d:%d\n\n",ii,aa);
			printf("Characteristic Length = %Lf\n A/M = %Lf\n Cross Sectional Area = %Lf\n",
				list[ii]->head[aa]->charLen,list[ii]->head[aa]->areaMass,
					list[ii]->head[aa]->csArea);
			printf("Mass = %Lf\n Delta V = %Lf\n",list[ii]->head[aa]->mass,
				list[ii]->head[aa]->deltaV);
			printf("**************************************\n\n");
			aa++;
		}
		aa=0;
		ii++;
	}
}
void printFile(linkedList** list)
{
	int ii=0,aa=0,numNodes,numBins;
	char* fn = "fragment.dat";
	FILE *fp;
	numBins = list[ii]->numbBins;
	fp = xfopen(fn,"w");
	while(ii<numBins)
	{
		numNodes=list[ii]->numbNodes;
		while(aa<numNodes)
		{
			fprintf(fp,"%Lf %Lf %Lf %Lf %Lf\n",list[ii]->head[aa]->charLen,
				list[ii]->head[aa]->areaMass,list[ii]->head[aa]->csArea,
				list[ii]->head[aa]->mass,list[ii]->head[aa]->deltaV);
			aa++;
		}
		aa=0;
		ii++;
	}
	xfclose(fp);
}
long double charLenUp(long double charLLen,double mass)
{
	long double fragBin,newUpper=charLLen;
	do
	{
		fragBin = 0.1*pow(mass,0.75)*(pow(newUpper,-1.71)-pow(newUpper+0.01,-1.71));
		newUpper+=0.01;
	}while(fragBin>1);	
	newUpper-=0.01;
	return newUpper;
}	
long double sampleCDF(long double** cumulDF)
{
	int ii =0, notFound = 1; 
	long double val,random;
	random=(1.0*rand()/2147483647.0);
	//printf("random = %Lf\n",random);
	while(notFound)
	{
		if(random<=cumulDF[1][ii])
		{
			notFound=0;
			val=cumulDF[0][ii];
			//printf("sample[%d] = %Lf\n",ii, val);
		}
		ii++;
	}
	return pow(10,val);
}
void createCDF(long double** probDF,long double** cumulDF, int xNum,long double lower,long double upper)
{
	long double iteration,cumulVal=0.0;
	int ii=4,aa=0;
	iteration = (upper-lower)/xNum;
	while(ii<=xNum-4)
	{
		cumulVal+=(2.0/45.0)*iteration*(7*probDF[1][ii-4] + 32*probDF[1][ii-3] +
			 12*probDF[1][ii-2] + 32*probDF[1][ii-1] + 7*probDF[1][ii]);
		cumulDF[1][aa]=cumulVal;
		aa++;
		ii+=4; 
	}
	cumulDF[1][aa]=cumulDF[1][aa-1];
	aa=0;
	ii=xNum/4;
	while(aa<ii)
	{
		cumulDF[1][aa]=cumulDF[1][aa]*(1/cumulDF[1][(xNum/4)-1]);
		aa++;
	}
}
void setStanVal(long double** plotXYval,int xNum,long double startX,long double stopX)
{
	long double iteration;
	int ii=0;
	iteration = (stopX-startX)/xNum;
	while(ii!=xNum)
	{
		plotXYval[0][ii]=startX;
		startX+=iteration;
		ii++;
	}
}
void addVelVal(long double** plotXYval,long double logAM,int xNum,int impactType)
{
	int ii=0;
	long double iteration,mu,xVal;
	if(impactType==0)
	{
		mu=0.2*logAM+1.85;
	}
	else
	{
		mu = 0.9*logAM+2.9;
	}
	iteration = 9.0/xNum;
	xVal=-2;
	while(ii!=xNum)
	{
		plotXYval[0][ii]=xVal;
		plotXYval[1][ii]=(1.0/(0.4*pow((2*M_PI),0.5)))*exp((-1*pow((plotXYval[0][ii]-mu),2)/
			(2*pow((0.4),2))));
		xVal+=iteration;
		ii++;
		//printf("check3.%d\n",ii);
	}
}
void addVal(long double** plotXYval,long double* parArray,int xNum)
{
	int ii=0;
	while(ii!=xNum)
	{
		plotXYval[1][ii]=parArray[0]*((1/(parArray[2]*pow((2*M_PI),0.5)))*
			exp((-1*pow((plotXYval[0][ii]-parArray[1]),2)/
			(2*pow((parArray[2]),2)))))+(1-parArray[0])*
			((1/(parArray[4]*pow((2*M_PI),0.5)))
			*exp((-1*pow((plotXYval[0][ii]-parArray[3]),2)
			/(2*pow((parArray[4]),2)))));
		ii++;
		//printf("check3.%d\n",ii);
	}
}
void normalise(long double** plotXYval,int xNum)
{
	long double largeVal;
	int ii=0;
	xNum--;
	largeVal=plotXYval[1][0];
	while(ii!=xNum)
	{
		if(largeVal<plotXYval[1][ii+1])
		{
			largeVal=plotXYval[1][ii+1];
		}
		ii++;
	}
	ii=0;
	xNum++;
	while(ii!=xNum)
	{
		plotXYval[1][ii]/=largeVal;
		ii++;
	}
}
void createDArray(int arSize,int xNum,long double ** array)
{
	int ii=0;
	while(ii!=arSize)
	{
		array[ii]=(long double*)calloc(xNum,sizeof(*array[ii]));
		ii++;
	}
}
void deAllocA(long double ** array,int length)
{
	int ii;
	for(ii=0;ii<length;ii++)
	{
		free(array[ii]);
	}
	free(array);
}
void deAllocLL(linkedList** list)
{
	int numBins,ii=0,aa=0;
	numBins=list[0]->numbBins;
	while(ii<numBins)
	{
		while(aa<list[ii]->numbNodes)
		{
			free(list[ii]->head[aa]);
			aa++;
		}
		aa=0;
		free(list[ii]->head);
		free(list[ii]);
		ii++;
	}
	free(list);
}
void createEq(long double* parArray)
{
	char str[200],strCom[200]="gnuplot -p -e 'plot  ";
	sprintf(str,"%Lf*((1/(%Lf*(2*pi)**0.5))*exp((-(x-%Lf)**2/(2*(%Lf)**2))))+(1-%Lf)*((1/(%Lf*(2*pi)**0.5))*exp((-(x-%Lf)**2/(2*(%Lf)**2))))",
			 parArray[0],parArray[2],parArray[1],parArray[2],parArray[0],parArray[4],parArray[3],parArray[4]);
	strcat(strCom,str);
	printf("\nEquation String : \n\n%s\n\n",str);
	strcat(strCom,"'");
	//system(strCom);

}
void getRBParam(long double charLen,long double * parArray)
{
	parArray[2]=0.55;
	parArray[3]=-0.9;
	//parArray[0]
	if(charLen<=-1.4)	
		{parArray[0]=1;}
	else if(-1.4<charLen && charLen<0)
		{parArray[0]=1-0.3571*(charLen+1.4);}
	else if(charLen>=0)
		{parArray[0]=0.5;}
	//parArray[1]
	if(charLen<=-0.5)
		{parArray[1]=-0.45;}
	else if(-0.5<charLen && charLen<0)
		{parArray[1]=-0.45-0.9*(charLen+0.5);}
	else if(charLen>=0)
		{parArray[1]=-0.9;}
	//parArray[4]
	if(charLen<=-1.0)
		{parArray[4]=0.28;}
	else if(-1.0<charLen && charLen<0.1)
		{parArray[4]=0.28-0.1636*(charLen+1);}
	else if(charLen>=0.1)
		{parArray[4]=0.1;}
}
void getSCParam(long double charLen,long double * parArray)
{
	//parArray[0]
	if(charLen<=-1.95)	
		{parArray[0]=0.0;}
	else if(-1.95<charLen && charLen<0.55)
		{parArray[0]=0.3+0.4*(charLen+1.2);}
	else if(charLen>=0.55)
		{parArray[0]=1.0;}
	//parArray[1]
	if(charLen<=-1.10)
		{parArray[1]=-0.6;}
	else if(-1.10<charLen && charLen<0.0)
		{parArray[1]=-0.6-0.318*(charLen+1.1);}
	else if(charLen>=0.0)
		{parArray[1]=-0.95;}
	if(charLen<=-1.30)
		{parArray[2]=0.1;}
	else if(-1.30<charLen && charLen<-0.30)
		{parArray[2]=0.1+0.2*(charLen+1.3);}
	else if(charLen>=-0.30)
		{parArray[2]=0.3;}
	//parArray[3]
	if(charLen<=-0.70)
		{parArray[3]=-1.2;}
	else if(-0.70<charLen && charLen<-0.10)
		{parArray[3]=-1.2-(4/3)*(charLen+0.7);}
	else if(charLen>=-0.10)
		{parArray[3]=-2;}
	//parArray[4]
	if(charLen<=-0.50)
		{parArray[4]=0.5;}
	else if(-0.50<charLen && charLen<-0.30)
		{parArray[4]=0.5-(charLen+0.5);}
	else if(charLen>=-0.30)
		{parArray[4]=0.3;}
	
}
int createScript(long double** plotXYval,int xNum,char* dataOne,char* dataTwo)
{
	int pid, status,ii;
	char *fn = "Model.dat";
	char fnbase[256] = "", fnplt[256] = "";
	FILE *fp = NULL;

	fp = xfopen (fn, "w");      /* open output file */

	for (ii = 0; ii < xNum; ii++)  /* write values to file */
	{
	        fprintf(fp, "%Lf %Lf\n", plotXYval[0][ii], plotXYval[1][ii]);
	}
	xfclose (fp);   /* close output file */

    /* create 'plot' file 'fn.plt' */
	strcpy (fnplt, fnwoext (fnbase, fn));
	strcat (fnplt, ".plt");
	if (!xfexists (fnplt))
	{
	        xfopen (fnplt, "w");
	        fprintf (fp, "set xlabel 'log(A/M)'\n"
                    "set ylabel 'Relative number of objects per A/M bin'font ',12'\n"
                    "set title 'A/M Distribution' font ',12'\n"
                    "set grid\n"
                    "set style data lines\n"
                    "plot \"%s\" using 1:2 lw 3 linecolor rgb \"blue\""
                    ,fn);
		if(dataOne && dataTwo)
		{
			fprintf(fp,", \\\n"
				"\"%s\" using 1:2 lw 1 linecolor rgb \"red\", \\\n"
				"\"%s\" using 1:2 lw 2 linecolor rgb \"green\"\n"
				"quit\n",dataOne,dataTwo);
		}
		else
		{
			fprintf(fp,"\nquit\n");
		}
	        xfclose (fp);
	}

    /* fill arguments array for execvp */
	char *args[] = { "gnuplot", "-p", fnplt, NULL };

	if ((pid = (fork())) < 0) 
	{ /* fork plot process */
		fprintf (stderr, "fork() error: fork failed.\n");
        	return 1;
    	}
	else if (pid == 0)
	{    /* plot from child process */
	        if (execvp (*args, args) == -1)
		{
	        	fprintf (stderr, "execvp() error: returned error.\n");
			_exit (EXIT_FAILURE);
        	}
	}

	waitpid (pid, &status, 0);  /* wait for plot completion (not req'd) */

	return 0;
}

/** fopen with error checking - short version */
FILE *xfopen (const char *fn, const char *mode)
{
    if (!fn || !mode || badmode (mode)) {
        fprintf (stderr, "xfopen() error: invalid parameter.\n");
        exit (EXIT_FAILURE);
    }
    FILE *fp = fopen (fn, mode);

    if (!fp) {
        fprintf (stderr, "xfopen() error: file open failed '%s'.\n", fn);
        exit (EXIT_FAILURE);
    }

    return fp;
}

/** validate file mode 's' is "rwa+b" */
int badmode (const char *s)
{
    const char *modes = "rwa+b";

    for (; *s; s++) {
        const char *m = modes;
        int valid = 0;
        while (*m) if (*s == *m++) { valid = 1; break; }
        if (!valid) return *s;
    }
    return 0;
}

/** file close with error check */
int xfclose (FILE *fp)
{
    if (fclose (fp)) {
        fprintf (stderr, "xfclose() error: nonzero return on fclose.\n");
        return 1;
    }
    return 0;
}

/** check if file 'fn' already exists */
int xfexists (char *fn)
{
    /* if access return is not -1 file exists */
    if (access (fn, F_OK ) != -1 )
        return 1;

    return 0;
}

/** isolate filename, without path or extension */
char *fnwoext (char *nm, char *fn)
{
    char *p  = NULL, *ep = NULL;
    char fnm[256] = "";

    if (!fn) return NULL;
    strcpy (fnm, fn);
    if ((p = strrchr (fnm, '/')))
        p++;
    else
        p = fnm;

    if ((ep = strrchr (p, '.'))) {
        *ep = 0;
        strcpy (nm, p);
        *ep = '.';
    } else
        strcpy (nm, p);

    return nm;
}	
			
		
