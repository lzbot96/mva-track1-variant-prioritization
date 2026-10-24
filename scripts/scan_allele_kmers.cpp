#include <zlib.h>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
struct Pattern {std::string name,seq; uint64_t code;};
int base(char c){switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}}
int main(int argc,char**argv){
 if(argc!=5){std::cerr<<"usage: patterns.tsv input.fastq.gz counts.tsv matched.fastq\n";return 2;}
 std::ifstream pf(argv[1]);std::vector<Pattern> pats;std::array<std::vector<int>,65536> lookup;std::string line;
 while(std::getline(pf,line)){std::istringstream ss(line);Pattern p;ss>>p.name>>p.seq;if(p.seq.size()!=31)return 3;p.code=0;for(char c:p.seq){int b=base(c);if(b<0)return 3;p.code=(p.code<<2)|b;}lookup[p.code&65535].push_back(pats.size());pats.push_back(p);}
 if(pats.empty())return 3;
 gzFile in=gzopen(argv[2],"rb");if(!in)return 4;gzbuffer(in,4*1024*1024);
 std::ofstream matches(argv[4]);if(!matches){gzclose(in);return 7;}std::vector<uint64_t> counts(pats.size(),0),q20(pats.size(),0),q30(pats.size(),0);uint64_t records=0;int row=0;std::string header,seq,pending;std::vector<std::pair<int,int>> hits;
 auto process=[&](const std::string &s){
  if(row==0){header=s;if(s.empty()||s[0]!='@')throw std::runtime_error("invalid header");}
  else if(row==1){seq=s;hits.clear();uint64_t code=0;int valid=0;for(size_t i=0;i<s.size();i++){int b=base(s[i]);if(b<0){valid=0;code=0;continue;}code=((code<<2)|b)&0x3fffffffffffffffULL;valid++;if(valid>=31){for(int j:lookup[code&65535])if(code==pats[j].code)hits.push_back({j,(int)i-15});}}}
  else if(row==2){if(s.empty()||s[0]!='+')throw std::runtime_error("invalid separator");}
  else {if(s.size()!=seq.size())throw std::runtime_error("sequence/quality length mismatch");records++;std::vector<bool> seen(pats.size(),false);for(auto h:hits){int j=h.first;if(!seen[j]){counts[j]++;int q=(unsigned char)s[h.second]-33;q20[j]+=q>=20;q30[j]+=q>=30;seen[j]=true;}}if(!hits.empty())matches<<header<<'\n'<<seq<<"\n+\n"<<s<<'\n';}
  row=(row+1)%4;
 };
 std::vector<char> buffer(4*1024*1024);int n;
 try{while((n=gzread(in,buffer.data(),buffer.size()))>0){size_t start=0;for(int i=0;i<n;i++)if(buffer[i]=='\n'){pending.append(buffer.data()+start,i-start);if(!pending.empty()&&pending.back()=='\r')pending.pop_back();process(pending);pending.clear();start=i+1;}pending.append(buffer.data()+start,n-start);}if(n<0)throw std::runtime_error("gzip read failed");if(!pending.empty())process(pending);if(row!=0)throw std::runtime_error("truncated FASTQ");}catch(const std::exception&e){std::cerr<<e.what()<<'\n';gzclose(in);return 5;}
 if(gzclose(in)!=Z_OK)return 6;
 matches.close();if(!matches)return 7;
 std::ofstream out(argv[3]);if(!out)return 7;out<<"pattern\treads\tq20_reads\tq30_reads\tfile_records\n";for(size_t i=0;i<pats.size();i++)out<<pats[i].name<<'\t'<<counts[i]<<'\t'<<q20[i]<<'\t'<<q30[i]<<'\t'<<records<<'\n';out.close();if(!out)return 7;std::cerr<<"Processed "<<records<<" records\n";
}
