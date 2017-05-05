// samtools mpileup NA06985.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam | ./pileup2shortform -outfile temp -fai index.fai
// samtools mpileup ../test/smallBam/smallNA07000.mapped.ILLUMINA.bwa.CEU.low_coverage.20111114.bam | ./pileup2shortform -outfile temp -fai /home/rowdy/github/angsd/test/hg19.fa.fai
// g++ -O3 -o pileup2shortform pileup2shortform.cpp

#include <iostream>

#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>
#include <sys/stat.h>
#include <cstring>
#include <string>
#include <map>
#include <utility>
#include <unistd.h>
#define LENS 10000000
#define STRLEN 1000 //max number of chars per outputfilename

//std::string oFilename = "outfile.txt";
const char *oFilename = NULL;
char *outFile = new char[STRLEN];
const char* qfilename = "argfile.txt";
const char* faifilename = "fai";

typedef struct{
  const char * chr;
  int len;
  //  int u1; the unknown stuff1
  //  int u2; the unknown stuff2
}aFai;
struct cmp_str
{
   bool operator()(char const *a, char const *b)
   {
      return std::strcmp(a, b) < 0;
   }
};

typedef std::map<const char*,aFai,cmp_str> aMap;


aMap getFai(const char *fname){
  FILE *fp = NULL;
  fp=fopen(fname,"r");
  if(fp==NULL){
    fprintf(stderr,"problems opening: %s\n",fname);
    exit(0);
  }
  aMap myMap;
  char buf[10000];
  while(fgets(buf,10000,fp)){
    char *chr = strtok(buf,"\n \t");
    aFai ret;
    ret.chr = strdup(chr);
    ret.len = atoi(strtok(NULL,"\n \t"));
    myMap.insert(std::pair<char*,aFai>(strdup(chr),ret));
  }
  fprintf(stderr,"number of chr in fai:%s = %lu\n",fname,myMap.size());
  return myMap;
}


int main(int argc, char *argv[]){
  int useQ=0;
  int argPos = 1;
  int maxQ=1000;
  int minQ=0;
  int offset = 33;
  aMap myMap;
  while(argPos <argc){
    if (strcmp(argv[argPos],"-outfile")==0){
      oFilename = argv[argPos+1];   
    }
    else if(strcmp(argv[argPos],"-qfile")==0){
      useQ=1;
      fprintf(stdout,"using custom quality scores\n");
      qfilename = argv[argPos+1]; 
    }
    else if(strcmp(argv[argPos],"-maxQ")==0){
      maxQ= atoi(argv[argPos+1]);
    }
    else if(strcmp(argv[argPos],"-fai")==0){
      faifilename = argv[argPos+1]; 
      myMap = getFai(faifilename);
    }
    else if(strcmp(argv[argPos],"-minQ")==0){
      minQ= atoi(argv[argPos+1]);
    }
    else if(strcmp(argv[argPos],"-off")==0){
     offset = atoi(argv[argPos+1]); 
    }
 
    //    if (strcmp(argv[argPos],"-qfile")==0){
    //      oFilename = argv[argPos+1];   
    //  }
    argPos+=2;
  }
  if(isatty(fileno(stdin))){
    fprintf(stderr,"You should really pipe data into the program\n");
    return 0;
  }
  
 

  const char *delim = " \n\t";
  char *buf =(char *) malloc(LENS);
  char *chr;
  int pos;
  char *seq;
  char *qual;
  

  int qCutoff[9];
  for(int i=0;i<9;i++)
    qCutoff[i] = minQ;
 

  if(useQ){
    fprintf(stdout,"quality cutOff \n");
    FILE *infileP=NULL;
    infileP=fopen(qfilename,"r");
    int line=0;
    while(fgets(buf,LENS,infileP)){
      qCutoff[line]=atoi(strtok(buf,delim));
      fprintf(stdout,"%d\n",qCutoff[line]);
      fflush(stdout);
      line++;
    }
    fclose(infileP);
  }
  int lookup[255];
  for(int i=0;i<255;i++)
    lookup[i]=-1;
  lookup['A'] = 0;
  lookup['C'] = 1;
  lookup['G'] = 2;
  lookup['T'] = 3;
  lookup['N'] = 4;
  lookup['a'] = 5;
  lookup['c'] = 6;
  lookup['g'] = 7;
  lookup['t'] = 8;
  lookup['n'] = 4;

  char intToRef[9] = {'A','C','G','T','N','a','c','g','t'};

  char *current_chr=new char[STRLEN];
  FILE *outfileP=NULL;

  size_t nLines=0;

  int genomePos=1;
  while(fgets(buf,LENS,stdin)){
    //    nLines++;
    //if(nLines>10000)
    //  break;
    int count[4]={0,0,0,0};
    chr = strtok(buf,delim);

    if(outfileP==NULL){ //first chromosome
      fprintf(stderr,"starting chromsome: %s\n",chr);
      aMap::iterator it = myMap.find(chr);
      if(it==myMap.end()){
	fprintf(stderr,"chr: %s doesn't exist in faifail:%s\n",chr,faifilename);
	exit(0);
      }
      aFai fai = it->second;
      //      fai.len;
      fprintf(stderr,"starting chromsome: %s with refLen: %d\n",chr,fai.len);
      snprintf(outFile,STRLEN,"%s.%s",oFilename,chr);
      outfileP = fopen(outFile,"w");
      strcpy(current_chr,chr);
    }
    //    fprintf(stderr,"val %d %s %s\n",strcmp(chr,current_chr),chr,current_chr);
    //    fprintf(stderr,"strlen1=%d\tstrlen2=%d\n",strlen(chr),strlen(current_chr));
    if(strcmp(chr,current_chr)!=0){ // new chromosome
      //print rest of chromosome
      aMap::iterator oldMap = myMap.find(current_chr);
      aFai oldfai = oldMap->second;
      fprintf(stderr,"printing rest of chr: %s with refLen: %d, genomePos: %d\n",current_chr,oldfai.len,genomePos);
      while(genomePos<=oldfai.len){
      //      fprintf(outfileP,"%d\tN\n",genomePos);
	fprintf(outfileP,"N");
	genomePos++;
      }


      aMap::iterator it = myMap.find(chr);
      if(it==myMap.end()){
	fprintf(stderr,"chr: %s doesn't exist in faifail:%s\n",chr,faifilename);
	exit(0);
      }
      aFai fai = it->second;
      //      fai.len;
      fprintf(stderr,"starting chromsome: %s with refLen: %d\n",chr,fai.len);
      fflush(stderr);
      if(outfileP!=NULL)
	fclose(outfileP);
      strcpy(current_chr,chr);
      snprintf(outFile,STRLEN,"%s.%s",oFilename,chr);
      outfileP = fopen(outFile,"w");
      genomePos=1;
    }

    pos = atoi(strtok(NULL,delim));
    char ref = strtok(NULL,delim)[0];
    //strip 2 columns
    strtok(NULL,delim);
    seq = strtok(NULL,delim);
    qual = strtok(NULL,delim);


    int base=4;
    uint j=0;
    for(uint i=0;i<strlen(qual);i++){
      while(seq[j]=='^'||seq[j]=='$'){
	if(seq[j]=='^'){
	  j++;
	  j++;
	}
	else if(seq[j]=='$'){
	  j++;
	}
      }
 
      if(seq[j]=='+'||seq[j]=='-'||seq[j]=='*')
	break;
      int tmp=lookup[seq[j]];
      if(qual[i]<(offset+qCutoff[tmp])||qual[i]>(offset+maxQ)){
	j++;
	continue;
      }
      if(seq[j]=='.'||seq[j]==',')
	tmp=lookup[ref];

      if(tmp==-1){
	fprintf(stderr,"error: pos %d chr %s is not valid %c\n",pos,chr,seq[j]);
	fflush(stderr);
	exit(0);
      }
      if(tmp!=4){
	count[tmp%5]++;
      }

      if(tmp!=4&&base==4){
	base=tmp;

      }
      j++;
    } //for
    if(seq[j]=='$')
      j++;
    else if(seq[j]=='^')
      j=j+2;
    if(j-strlen(seq)!=0){
      if(seq[j]=='+'||seq[j]=='-'||seq[j]=='*')
	continue;
      fprintf(stderr,"error: pos %d chr %s lengthFUck %d\n",pos,chr,j);
    }


    while(genomePos<pos){
      //      fprintf(outfileP,"%d\tN\n",genomePos);
      fprintf(outfileP,"N");
      genomePos++;
    }

    fprintf(outfileP,"%c",intToRef[base]);
    genomePos++;


    // fwrite(&pos,sizeof(int),1,outfileP);
    //  fwrite(count,sizeof(int)*4,1,outfileP);
    //    fprintf(outfileP,"%d\n",pos);
    
  } //while


  //print out rest of last chromosome
  aMap::iterator oldMap = myMap.find(current_chr);
  aFai oldfai = oldMap->second;
  while(genomePos<=oldfai.len){
    //      fprintf(outfileP,"%d\tN\n",genomePos);
    fprintf(outfileP,"N");
    genomePos++;
  }


  return 0;
}
