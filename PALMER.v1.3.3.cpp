 //copyright by ArthurZhou @ UMich&FUDAN&HUST
#include <stdlib.h>
#include <iostream>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <cstdlib>
#include <numeric>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "scp/tube.cpp"

using namespace std;

int main(int argc, char *argv[]){

//parameters_start
    
    string T, WD, inputF, output, SP, ref, CHR, ref_fa, cus, cusin;
    //int NUM_threads=30;
    int NUM_threads=10;
    ifstream file1;
    ifstream file11;
    //ifstream file12;
    ifstream file13;
    //ofstream file2;
    int flag_wd=0;
    int flag_inputf=0;
    int flag_T=0;
    int flag_reffa=0;
    int flag_cus=0;
    int flag_cusin=0;
    output="output.txt";
    SP="Human";
    //T="LINE";
    CHR="ALL";
    int ref_n=0;
    
    string dir;
    dir=argv[0];
    
    for(int i=1;i!=argc;i++){
        if(strncmp(argv[i],"--workdir",9)==0){
            WD=argv[i+1];
            flag_wd=1;
        }
        if(strncmp(argv[i],"--chr",5)==0){
            CHR=argv[i+1];
        }
        if(strncmp(argv[i],"--ref_ver",9)==0){
            ref=argv[i+1];
            if(ref=="GRCh37") ref_n=37;
            else if(ref=="GRCh38") ref_n=38;
            else if(ref=="hg19") ref_n=19;
            else {
                cout<<"PLEASE INPUT A CORRECT HUMAN REFERENCE :("<<endl;
            }
        }
        if(strncmp(argv[i],"--type",6)==0){
            T=argv[i+1];
            flag_T=1;
        }
        if(strncmp(argv[i],"--input",7)==0){
            flag_inputf=1;
            file1.open(argv[i+1],ios::in | ios::binary);
            if (!file1.is_open())
            {
                cout <<"CANNOT OPEN INPUT FILE :("<< endl;
                continue;
            }
            
            inputF=argv[i+1];
        }
        if(strncmp(argv[i],"--output",7)==0){
            output=argv[i+1];
        }
        if(strncmp(argv[i],"--species",9)==0){
            SP=argv[i+1];
            cout <<"WE DO NOT NEED SPECIES PARAMETER RIGHT NOW :)"<< endl;
        }
        if(strncmp(argv[i],"--thread",8)==0){
            //NUM_threads=argv[i+1];
            cout <<"WE DO NOT SUPPORT MULTITHREADS RIGHT NOW :)"<< endl;
        }
        if(strncmp(argv[i],"--ref_fa",8)==0){
            flag_reffa=1;
            file11.open(argv[i+1],ios::in | ios::binary);
            if (!file11.is_open())
            {
                cout <<"CANNOT OPEN REFERENCE FILE :("<< endl;
                cout<<"PLEASE ASSIGN A REFERENCE FILE."<<endl;
                //exit(1);
                flag_reffa=0;
            }
            ref_fa=argv[i+1];
        }
        
        if(strncmp(argv[i],"--custom_seq",12)==0){
            flag_cus=1;
            /*
            file12.open(argv[i+1],ios::in | ios::binary);
            if (!file12.is_open())
            {
                cout <<"CANNOT OPEN CUSTOMIZED FILE :("<< endl;
                cout<<"PLEASE ASSIGN A CUSTOMIZED FASTA FILE."<<endl;
                //exit(1);
                flag_cus=0;
            }*/
            cus=argv[i+1];
            T="CUSTOMIZED";
        }
        if(strncmp(argv[i],"--custom_index",14)==0){
            flag_cusin=1;
            if(flag_cus==1){
                 file13.open(argv[i+1],ios::in | ios::binary);
                 if (!file13.is_open())
                 {
                 cout <<"You are not assigning any index files"<< endl;
                 cout<<"Custmized finding will initiate without masking module."<<endl;
                 //exit(1);
                 flag_cusin=0;
                 }
            }
            cusin=argv[i+1];
        }
        
    }
    
    
    if((flag_T==0&&flag_cus==0)||(flag_T==1&&flag_cus==1)||flag_wd==0||flag_inputf==0||ref_n==0||flag_reffa==0){
        if(flag_T==0&&flag_cus==0){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE! LINE/ALU/SVA"<<endl;}
        if(flag_T==1&&flag_cus==1){
            cout<<"***ERROR*** PLEASE ASSIGN A MEI TYPE WITHOUT YOUR CUSTOMIZED SEQUENCE"<<endl;}
        if(flag_wd==0){
            cout<<"***ERROR*** PLEASE SET UP A WORKING DIRECOTY!"<<endl;}
        if(flag_inputf==0){
            cout<<"***ERROR*** PLEASE INPUT A FILE!"<<endl;}
        if(ref_n==0||flag_reffa==0){
            cout<<"***ERROR*** PLEASE ASSIGN A CORRECT REFERENCE version/fasta!"<<endl;}
        cout<<endl;
        cout<<"***PALMER:Pre-mAsking Long reads for Mobile Element inseRtion***"<<endl;
        cout<<"Version: Beta1.2"<<endl;
        cout<<"Presented by Weichen Zhou @ Mills Lab. Sep.5th.2018"<<endl;
        cout<<endl;
        cout<<"Usage:"<<endl;
        cout<<endl;
        cout<<"--input"<<endl;
        cout<<"         aligned long-read sequencing BAM file with directory path"<<endl;
        cout<<endl;
        cout<<"--workdir"<<endl;
        cout<<"         the user's working directory"<<endl;
        cout<<endl;
        cout<<"--ref_ver (options: hg19, GRCh37 or GRCh38)"<<endl;
        cout<<"         reference genome used for the aligned file (only human genome rightnow)"<<endl;
        cout<<endl;
        cout<<"--ref_fa"<<endl;
        cout<<"         indexed fasta file of reference genome fasta file with directory path used for the aligned file"<<endl;
        cout<<endl;
        cout<<"--type (options: LINE, ALU, SVA, or coustomized (if you want to setup your costomized sequence))"<<endl;
        cout<<"         type of MEIs to detect"<<endl;
        cout<<endl;
        cout<<"--chr (default: whole genome; options: chr1, chr2, ...chrY)"<<endl;
        cout<<"         chr name for PALMER to run (if running for whole genome, don't need to assign)"<<endl;
        cout<<endl;
        
        cout<<"--custom_seq (default:no input)"<<endl;
        cout<<"         fasta file with directory path to customize your insertion finding"<<endl;
        cout<<endl;
        cout<<"--custom_index (default:no input, if you have --custom_seq parameter without --custom_index, PALMER will work without masking step)"<<endl;
        cout<<"         index file with directory path to mask the genome for your insertion finding"<<endl;
        cout<<endl;
        
        
        cout<<"--output (default: output)"<<endl;
        cout<<"         prefix of output file"<<endl;
        cout<<endl;
        exit(1);
    }
    
    
    string parameter[6];
    parameter[0]=T;
    parameter[1]=WD;
    parameter[2]=inputF;
    parameter[3]=output;
    parameter[4]=CHR;
    parameter[5]=ref;
    
    cout<<"Variant type is "<<parameter[0]<<endl;
    cout<<"Working directory is "<<parameter[1]<<endl;
    cout<<"Input file is "<<parameter[2]<<endl;
    cout<<"Output file is "<<parameter[1]<<parameter[3]<<endl;
    cout<<"Running on "<<parameter[4]<<endl;
    cout<<"ref is "<<parameter[5]<<endl;
    
//Buildup & index or HG version chose
    string sys_dir="dirname "+dir;
    char *syst_dir=new char[sys_dir.length()+1];
    strcpy(syst_dir, sys_dir.c_str());
    
    vector<string> dir_conv;
    dir_conv.clear();
    FILE *pp =popen(syst_dir,"r");
    char tmp[1024];
    while (fgets(tmp, sizeof(tmp), pp) != NULL) {
        if (tmp[strlen(tmp) - 1] == '\n') {
            tmp[strlen(tmp) - 1] = '\0';
        }
        dir_conv.push_back(tmp);
    }
    pclose(pp);
    string direc;
    direc=accumulate(dir_conv.begin(),dir_conv.end(),direc);
    
    string buildup=direc+"/index/";
    
    string sys_region_index,sys_line_region;
    
    if(T=="CUSTOMIZED"){
        
        //sys_region_index=buildup+"region.split.index.GRCh37";
        if(flag_cusin==1){
            sys_line_region=cusin;
        }
        else {
            sys_line_region="NULL";
        }
        
        T=cus;
    }
    
    
    if(ref_n==37){
        
        sys_region_index=buildup+"region.split.index.GRCh37";
        
        //test
        //sys_region_index=buildup+"region.split.index.test.txt";
        //
        
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh37";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh37";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh37";
        }
    }
    else if(ref_n==38){
        sys_region_index=buildup+"region.split.index.GRCh38";
        
        //test
        //sys_region_index=buildup+"region.split.index.test.txt";
        //
        
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.GRCh38";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.GRCh38";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.GRCh38";
        }
    }
    else if(ref_n==19){
        sys_region_index=buildup+"region.split.index.hg19";
        
        //test
        //sys_region_index=buildup+"region.split.index.test.txt";
        //
        
        if(T=="LINE"){
            sys_line_region=buildup+"LINEs.regions.hg19";
        }
        if(T=="ALU"){
            sys_line_region=buildup+"Alu.regions.hg19";
        }
        if(T=="SVA"){
            sys_line_region=buildup+"SVA.regions.hg19";
        }
    }
    //cout<<sys_region_index<<" "<<sys_line_region<<endl;
 //original

    //string sys_region_index=buildup+"region.split.index";
    char *syst_region_index =new char[sys_region_index.length()+1];
    strcpy(syst_region_index, sys_region_index.c_str());
    
    ifstream file2;
    file2.open(syst_region_index);
    
    if (!file2.is_open())
    {
        cout <<"CANNOT OPEN INDEX FILE"<< endl;
        //cout<<"YEA"<<endl;
        exit(1);
    }
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    /*
    string sys_makeref = "makeblastdb -in "+WD+"buildup/L1.3.fasta  -dbtype nucl -parse_seqids";
    char *syst_makeref = new char[sys_makeref.length()+1];
    strcpy(syst_makeref, sys_makeref.c_str());
    system(syst_makeref);
     */
//parameters_end
    
//multiple threads
    //##########hard code##########upgrade section##########
    
    string input_index;
    int line_index=0;
    for(int i=0;!file2.eof();){
        //file2>>input_index;
        file2>>input_index;
        if(CHR=="ALL"){
            line_index=i;
            i++;
        }
        else if(CHR=="chr1"&&input_index=="1"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr2"&&input_index=="2"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr3"&&input_index=="3"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr4"&&input_index=="4"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr5"&&input_index=="5"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr6"&&input_index=="6"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr7"&&input_index=="7"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr8"&&input_index=="8"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr9"&&input_index=="9"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr10"&&input_index=="10"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr11"&&input_index=="11"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr12"&&input_index=="12"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr13"&&input_index=="13"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr14"&&input_index=="14"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr15"&&input_index=="15"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr16"&&input_index=="16"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr17"&&input_index=="17"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr18"&&input_index=="18"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr19"&&input_index=="19"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr20"&&input_index=="20"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr21"&&input_index=="21"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chr22"&&input_index=="22"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chrX"&&input_index=="X"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chrY"&&input_index=="Y"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR=="chrM"&&input_index=="MT"&&ref_n==37){
            line_index=i;
            i++;
        }
        else if(CHR==input_index&&(ref_n==38||ref_n==19)){
            line_index=i;
            i++;
        }
        file2>>input_index;
        file2>>input_index;
        
    }
    
    line_index=line_index+1;
    
    cout<<"THERE ARE "<<line_index<<" REGIONS TO COUNT."<<endl;
    cout<<"Pre-masking step & single read calling step is initiated."<<endl;
    
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    int NUM_circle;
    NUM_circle=(line_index/(NUM_threads+1))+1;

    //pid_t p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
    //, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26, p27, p28, p29, p30;
    
    string input;
    string chr;
    int start, end;

//no multithread
    for(int i=0;i!=line_index;){
        
        file2>>chr;
        file2>>start;
        file2>>end;
        int chr_index=0;
        if(CHR=="ALL"){
            chr_index=1;
        }
        else if(CHR=="chr1"&& chr=="1"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr2"&& chr=="2"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr3"&& chr=="3"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr4"&& chr=="4"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr5"&& chr=="5"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr6"&& chr=="6"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr7"&& chr=="7"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr8"&& chr=="8"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr9"&& chr=="9"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr10"&& chr=="10"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr11"&& chr=="11"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr12"&& chr=="12"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr13"&& chr=="13"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr14"&& chr=="14"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr15"&& chr=="15"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr16"&& chr=="16"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr17"&& chr=="17"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr18"&& chr=="18"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr19"&& chr=="19"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr20"&& chr=="20"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr21"&& chr=="21"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chr22"&& chr=="22"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chrX"&& chr=="X"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chrY"&& chr=="Y"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR=="chrM"&& chr=="MT"&&ref_n==37){
            chr_index=1;
        }
        else if(CHR==chr&&(ref_n==38||ref_n==19)){
            chr_index=1;
            //cout<<"right call"<<endl;
        }
        if(chr_index==1){
            i++;
            //cout<<chr<<endl;
            //getchar();
            //****
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_fa);
        }
        /*ver1.3
        if(chr_index==1){
            tube(WD, inputF, chr, start, end, sys_line_region, T, ref_n, direc, ref_file);
        }*/
    }
 
 
//merge and calling

    cout<<"Merging step is initiated."<<endl;
    //mkdir
    
    string WD_chr1=WD+"chr1/";
    string WD_chr2=WD+"chr2/";
    string WD_chr3=WD+"chr3/";
    string WD_chr4=WD+"chr4/";
    string WD_chr5=WD+"chr5/";
    string WD_chr6=WD+"chr6/";
    string WD_chr7=WD+"chr7/";
    string WD_chr8=WD+"chr8/";
    string WD_chr9=WD+"chr9/";
    string WD_chr10=WD+"chr10/";
    string WD_chr11=WD+"chr11/";
    string WD_chr12=WD+"chr12/";
    string WD_chr13=WD+"chr13/";
    string WD_chr14=WD+"chr14/";
    string WD_chr15=WD+"chr15/";
    string WD_chr16=WD+"chr16/";
    string WD_chr17=WD+"chr17/";
    string WD_chr18=WD+"chr18/";
    string WD_chr19=WD+"chr19/";
    string WD_chr20=WD+"chr20/";
    string WD_chr21=WD+"chr21/";
    string WD_chr22=WD+"chr22/";
    string WD_chr23=WD+"chrX/";
    string WD_chr24=WD+"chrY/";
    string WD_chr25=WD+"chrM/";
    
    string sys_WD_chr1="mkdir "+WD+"chr1/";
    string sys_WD_chr2="mkdir "+WD+"chr2/";
    string sys_WD_chr3="mkdir "+WD+"chr3/";
    string sys_WD_chr4="mkdir "+WD+"chr4/";
    string sys_WD_chr5="mkdir "+WD+"chr5/";
    string sys_WD_chr6="mkdir "+WD+"chr6/";
    string sys_WD_chr7="mkdir "+WD+"chr7/";
    string sys_WD_chr8="mkdir "+WD+"chr8/";
    string sys_WD_chr9="mkdir "+WD+"chr9/";
    string sys_WD_chr10="mkdir "+WD+"chr10/";
    string sys_WD_chr11="mkdir "+WD+"chr11/";
    string sys_WD_chr12="mkdir "+WD+"chr12/";
    string sys_WD_chr13="mkdir "+WD+"chr13/";
    string sys_WD_chr14="mkdir "+WD+"chr14/";
    string sys_WD_chr15="mkdir "+WD+"chr15/";
    string sys_WD_chr16="mkdir "+WD+"chr16/";
    string sys_WD_chr17="mkdir "+WD+"chr17/";
    string sys_WD_chr18="mkdir "+WD+"chr18/";
    string sys_WD_chr19="mkdir "+WD+"chr19/";
    string sys_WD_chr20="mkdir "+WD+"chr20/";
    string sys_WD_chr21="mkdir "+WD+"chr21/";
    string sys_WD_chr22="mkdir "+WD+"chr22/";
    string sys_WD_chr23="mkdir "+WD+"chrX/";
    string sys_WD_chr24="mkdir "+WD+"chrY/";
    string sys_WD_chr25="mkdir "+WD+"chrM/";
    
    if(CHR=="chr1"||CHR=="ALL"){
        char *syst_WD_chr1 =new char[sys_WD_chr1.length()+1];
        strcpy(syst_WD_chr1, sys_WD_chr1.c_str());
        system(syst_WD_chr1);
    }
    
    if(CHR=="chr2"||CHR=="ALL"){
        char *syst_WD_chr2 =new char[sys_WD_chr2.length()+1];
        strcpy(syst_WD_chr2, sys_WD_chr2.c_str());
        system(syst_WD_chr2);
    }
    
    if(CHR=="chr3"||CHR=="ALL"){
        char *syst_WD_chr3 =new char[sys_WD_chr3.length()+1];
        strcpy(syst_WD_chr3, sys_WD_chr3.c_str());
        system(syst_WD_chr3);
    }
    
    if(CHR=="chr4"||CHR=="ALL"){
        char *syst_WD_chr4 =new char[sys_WD_chr4.length()+1];
        strcpy(syst_WD_chr4, sys_WD_chr4.c_str());
        system(syst_WD_chr4);
    }
    
    if(CHR=="chr5"||CHR=="ALL"){
        char *syst_WD_chr5 =new char[sys_WD_chr5.length()+1];
        strcpy(syst_WD_chr5, sys_WD_chr5.c_str());
        system(syst_WD_chr5);
    }
    
    if(CHR=="chr6"||CHR=="ALL"){
        char *syst_WD_chr6 =new char[sys_WD_chr6.length()+1];
        strcpy(syst_WD_chr6, sys_WD_chr6.c_str());
        system(syst_WD_chr6);
    }
    
    if(CHR=="chr7"||CHR=="ALL"){
        char *syst_WD_chr7 =new char[sys_WD_chr7.length()+1];
        strcpy(syst_WD_chr7, sys_WD_chr7.c_str());
        system(syst_WD_chr7);
    }
    
    if(CHR=="chr8"||CHR=="ALL"){
        char *syst_WD_chr8 =new char[sys_WD_chr8.length()+1];
        strcpy(syst_WD_chr8, sys_WD_chr8.c_str());
        system(syst_WD_chr8);
    }
    
    if(CHR=="chr9"||CHR=="ALL"){
        char *syst_WD_chr9 =new char[sys_WD_chr9.length()+1];
        strcpy(syst_WD_chr9, sys_WD_chr9.c_str());
        system(syst_WD_chr9);
    }
    
    if(CHR=="chr10"||CHR=="ALL"){
        char *syst_WD_chr10 =new char[sys_WD_chr10.length()+1];
        strcpy(syst_WD_chr10, sys_WD_chr10.c_str());
        system(syst_WD_chr10);
    }
    
    if(CHR=="chr11"||CHR=="ALL"){
        char *syst_WD_chr11 =new char[sys_WD_chr11.length()+1];
        strcpy(syst_WD_chr11, sys_WD_chr11.c_str());
        system(syst_WD_chr11);
    }
    
    if(CHR=="chr12"||CHR=="ALL"){
        char *syst_WD_chr12 =new char[sys_WD_chr12.length()+1];
        strcpy(syst_WD_chr12, sys_WD_chr12.c_str());
        system(syst_WD_chr12);
    }
    
    if(CHR=="chr13"||CHR=="ALL"){
        char *syst_WD_chr13 =new char[sys_WD_chr13.length()+1];
        strcpy(syst_WD_chr13, sys_WD_chr13.c_str());
        system(syst_WD_chr13);
    }
    
    if(CHR=="chr14"||CHR=="ALL"){
        char *syst_WD_chr14 =new char[sys_WD_chr1.length()+1];
        strcpy(syst_WD_chr14, sys_WD_chr14.c_str());
        system(syst_WD_chr14);
    }
    
    if(CHR=="chr15"||CHR=="ALL"){
        char *syst_WD_chr15 =new char[sys_WD_chr15.length()+1];
        strcpy(syst_WD_chr15, sys_WD_chr15.c_str());
        system(syst_WD_chr15);
    }
    
    if(CHR=="chr16"||CHR=="ALL"){
        char *syst_WD_chr16 =new char[sys_WD_chr16.length()+1];
        strcpy(syst_WD_chr16, sys_WD_chr16.c_str());
        system(syst_WD_chr16);
    }
    
    if(CHR=="chr17"||CHR=="ALL"){
        char *syst_WD_chr17 =new char[sys_WD_chr17.length()+1];
        strcpy(syst_WD_chr17, sys_WD_chr17.c_str());
        system(syst_WD_chr17);
    }
    
    if(CHR=="chr18"||CHR=="ALL"){
        char *syst_WD_chr18 =new char[sys_WD_chr18.length()+1];
        strcpy(syst_WD_chr18, sys_WD_chr18.c_str());
        system(syst_WD_chr18);
    }
    
    if(CHR=="chr19"||CHR=="ALL"){
        char *syst_WD_chr19 =new char[sys_WD_chr19.length()+1];
        strcpy(syst_WD_chr19, sys_WD_chr19.c_str());
        system(syst_WD_chr19);
    }
    
    if(CHR=="chr20"||CHR=="ALL"){
        char *syst_WD_chr20 =new char[sys_WD_chr20.length()+1];
        strcpy(syst_WD_chr20, sys_WD_chr20.c_str());
        system(syst_WD_chr20);
    }
    
    if(CHR=="chr21"||CHR=="ALL"){
        char *syst_WD_chr21 =new char[sys_WD_chr21.length()+1];
        strcpy(syst_WD_chr21, sys_WD_chr21.c_str());
        system(syst_WD_chr21);
    }
    
    if(CHR=="chr22"||CHR=="ALL"){
        char *syst_WD_chr22 =new char[sys_WD_chr22.length()+1];
        strcpy(syst_WD_chr22, sys_WD_chr22.c_str());
        system(syst_WD_chr22);
    }
    
    if(CHR=="chrX"||CHR=="ALL"){
        char *syst_WD_chr23 =new char[sys_WD_chr23.length()+1];
        strcpy(syst_WD_chr23, sys_WD_chr23.c_str());
        system(syst_WD_chr23);
    }
    
    if(CHR=="chrY"||CHR=="ALL"){
        char *syst_WD_chr24 =new char[sys_WD_chr24.length()+1];
        strcpy(syst_WD_chr24, sys_WD_chr24.c_str());
        system(syst_WD_chr24);
    }
    
    if(CHR=="chrM"||CHR=="ALL"){
        char *syst_WD_chr25 =new char[sys_WD_chr25.length()+1];
        strcpy(syst_WD_chr25, sys_WD_chr25.c_str());
        system(syst_WD_chr25);
    }
   
    //##########hard code##########upgrade section##########end
    
    file2.close();
    file2.clear();
    file2.open(syst_region_index);
    
    
    //merge
    
    string sys_final_title = WD+output+"_calls.txt";
    char *syst_final_title = new char[sys_final_title.length()+1];
    strcpy(syst_final_title, sys_final_title.c_str());
    ofstream file3;
    file3.open(syst_final_title,ios::trunc);
    
    file3<<"cluster_id"<<'\t'<<"chr start1"<<'\t'<<"start2"<<'\t'<<"end1"<<'\t'<<"end2"<<'\t'<<"start1_inVariant"<<'\t'<<"start2_inVariant"<<'\t'<<"end1_inVariant"<<'\t'<<"end2_inVariant"<<'\t'<<"Confident_supporting_reads"<<'\t'<<"Potential_supporting_reads"<<'\t'<<"Ptential_segmental_supporting_reads"<<'\t'<<"orientation"<<'\t'<<"polyA-tail_size"<<'\t'<<"5'_TSD_size"<<'\t'<<"3'_TSD_size"<<'\t'<<"Predicted_transD_size"<<'\t'<<"Has_5'_inverted_sequence?"<<'\t'<<"5'_inverted_seq_end"<<'\t'<<"5'_seq_start"<<endl;
    
    string sys_final_tsd_title = WD+output+"_TSD_reads.txt";
    char *syst_final_tsd_title = new char[sys_final_tsd_title.length()+1];
    strcpy(syst_final_tsd_title, sys_final_tsd_title.c_str());
    ofstream file31;
    file31.open(syst_final_tsd_title,ios::trunc);
    
    file31<<"cluster_id"<<'\t'<<"read_name.info"<<'\t'<<"5'_TSD"<<'\t'<<"3'_TSD"<<'\t'<<"Predicted_transD"<<'\t'<<"Unique_26mer_at_5'junction"<<'\t'<<"Whole_insertion_seq"<<endl;
    
    for(int i=0;i!=line_index;){
        file2>>chr;
        file2>>start;
        file2>>end;
        
        stringstream ss1, ss2;
        ss1 << start;
        string s_start =ss1.str();
        ss2 << end;
        string s_end =ss2.str();
        
        string CHR_in;
        CHR_in="chr"+chr;

        if(ref_n==37&&CHR_in==CHR){
            
            string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
            //cout<<sys_final<<endl;
            char *syst_final = new char[sys_final.length()+1];
            strcpy(syst_final, sys_final.c_str());
            system(syst_final);
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            //cout<<sys_final<<endl;
            char *syst_final_tsd = new char[sys_final_tsd.length()+1];
            strcpy(syst_final_tsd, sys_final_tsd.c_str());
            system(syst_final_tsd);
            
            i++;
        }
        else if((CHR==chr)&&(ref_n==19||ref_n==38)){
           
            string sys_final="cat "+WD+chr+"_"+s_start+"_"+s_end+"/calls.txt >> "+sys_final_title;
            //cout<<sys_final<<endl;
            char *syst_final = new char[sys_final.length()+1];
            strcpy(syst_final, sys_final.c_str());
            system(syst_final);
            
            string sys_final_tsd="cat "+WD+chr+"_"+s_start+"_"+s_end+"/TSD_output.txt >> "+sys_final_tsd_title;
            //cout<<sys_final<<endl;
            char *syst_final_tsd = new char[sys_final_tsd.length()+1];
            strcpy(syst_final_tsd, sys_final_tsd.c_str());
            system(syst_final_tsd);
            
            i++;
        }
    }
    
    cout<<"Merging step completed."<<endl;
    cout<<"Final calls finished."<<endl;
    cout<<"Results are in "+WD+output<<endl;
    
}