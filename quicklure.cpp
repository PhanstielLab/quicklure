#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <tuple>
#include <unistd.h>
#include <unordered_set>
#include <vector>

const std::regex FASTA_REGEX(">?(chr(.+)):(\\d+)-(\\d+)");

void ScanProbesPass1(std::string candidatesFastaFilename, int repeatThresh, int GCThresh1, int GCThresh2, int GCThresh3, int GCThresh4,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueSitesFilename, std::string side,
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> &probes)
{
  std::ifstream fasta(candidatesFastaFilename);
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream is_file(issueSitesFilename, std::ofstream::app);
  if(fasta.good() && pl_file.good() && ip_file.good() && is_file.good())
  {
    std::string line;
    std::string description;
    int entries = 0;
    std::string sequence;
    std::string first_description;
    std::tuple<std::string, int, int, std::string, int, int, int> backup_probe;
    bool done = false;
    while(true)
    {
      char c = fasta.peek();
      if(c == '>' || c == EOF)
      {
        if(!description.empty())
        {
          if(entries == 1)
          {
            std::ofstream c1_file("condition1.txt", std::ofstream::app);
            if(c1_file.good())
              c1_file << description << ' ' << 0 << ' ' << 0 << ' ' << description << std::endl;
            c1_file.close();
            first_description = description;
          }
          int counter = 0;
          int C = 0;
          int G = 0;
          for(int i = 0; i < sequence.length(); i++)
          {
            if(islower(sequence[i]))
              counter++;
            if(sequence[i] == 'C')
              C++;
            if(sequence[i] == 'G')
              G++;
          }
          int GC = G + C;
          if(counter<=repeatThresh&&GC>=GCThresh1&&GC<=GCThresh2)
          {
            std::string chromosome;
            int start;
            int end;
            std::smatch matches;
            if(std::regex_search(description, matches, FASTA_REGEX))
            {
              chromosome = matches[1].str();
              start = stoi(matches[3].str());
              end = stoi(matches[4].str());
            }
            if(GC>=GCThresh3&&GC<=GCThresh4)
            {
              std::ofstream s("something.txt", std::ofstream::app);
              if(s.good())
                s << sequence << std::endl;
              s.close();
              probes.push_back(std::make_tuple(chromosome, start, end, sequence, sequence.length(), counter, GC));
              pl_file << description << ' ' << sequence << ' ' << sequence.length() << ' ' << counter << ' ' << GC << std::endl;
              done = true;
              break;
            }
            else if(std::get<0>(backup_probe).empty())
              backup_probe = std::make_tuple(chromosome, start, end, sequence, sequence.length(), counter, GC);
          }
          else
            ip_file << description << ' ' << counter << ' ' << GC << ' ' << sequence << std::endl;
        }
        if(c == EOF)
          break;
        description.clear();
        sequence.clear();
      }
      getline(fasta, line);
      if(c == '>')
      {
        description = line;
        entries++;
      }
      else
        sequence += line;
    }
    if(!done)
    {
      if(!std::get<0>(backup_probe).empty())
      {
        probes.push_back(backup_probe);
        pl_file << '>' << std::get<0>(backup_probe) << ':' << std::get<1>(backup_probe) << '-' << std::get<2>(backup_probe) << ' ' <<
          std::get<3>(backup_probe) << ' ' << std::get<4>(backup_probe) << ' ' << std::get<5>(backup_probe) << ' ' << std::get<6>(backup_probe) << std::endl;
      }
      else
      {
        std::smatch matches;
        if(std::regex_search(first_description, matches, FASTA_REGEX))
        {
          std::string chromosome = matches[1].str();
          int start = stoi(matches[3].str());
          int end = stoi(matches[4].str());
          is_file << "No " << side << " probes found for:" << chromosome << ' ' << (side == "upstream" ? end : start) << std::endl;
        }
      }
    }
    fasta.close();
    pl_file.close();
    ip_file.close();
    is_file.close();
  }
}

bool compareChromThenStart(const std::tuple<std::string, int, int, std::string, int, int, int> &left,
  const std::tuple<std::string, int, int, std::string, int, int, int> &right)
{
  return std::tie(std::get<0>(left), std::get<1>(left)) < std::tie(std::get<0>(right), std::get<1>(right));
}

void ScanProbesPass2_3(std::string candidatesFastaFilename, int repeatThresh, int GCThresh1, int GCThresh2,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueGapsFilename, int lineCounter,
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> &probes)
{
  std::ifstream fasta(candidatesFastaFilename);
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream ig_file(issueGapsFilename, std::ofstream::app);
  if(fasta.good() && pl_file.good() && ip_file.good() && ig_file.good())
  {
    std::string line;
    std::string description;
    std::string sequence;
    std::string chromosome;
    int pos1;
    int pos2;
    bool done = false;
    int end = 0;
    while(true)
    {
      char c = fasta.peek();
      if(c == '>' || c == EOF)
      {
        if(!description.empty())
        {
          if(pos1 > end)
          {
            int counter = 0;
            int C = 0;
            int G = 0;
            for(int i = 0; i < sequence.length(); i++)
            {
              if(islower(sequence[i]))
                counter++;
              if(sequence[i] == 'C')
                C++;
              if(sequence[i] == 'G')
                G++;
            }
            int GC = G + C;
            if (counter <= repeatThresh && GC >= GCThresh1 && GC <= GCThresh2)
            {
              probes.push_back(std::make_tuple(chromosome, pos1, pos2, sequence, sequence.length(), counter, GC));
              pl_file << chromosome << ' ' << pos1 << ' ' << pos2 << ' ' << sequence << ' ' << sequence.length() << ' ' << counter << ' ' << GC << std::endl;
              end = pos2;
              done = true;
            }
            else
              ip_file << description << ' ' << counter << ' ' << GC << std::endl;
            std::smatch matches;
            if(std::regex_search(description, matches, FASTA_REGEX))
            {
              chromosome = matches[1].str();
              pos1 = stoi(matches[3].str());
              pos2 = stoi(matches[4].str());
            }
          }
        }
        if(c == EOF)
          break;
        description.clear();
        sequence.clear();
      }
      getline(fasta, line);
      if(c == '>')
      {
        description = line;
        if(chromosome.empty())
        {
          std::smatch matches;
          if(std::regex_search(description, matches, FASTA_REGEX))
          {
            chromosome = matches[1].str();
            pos1 = stoi(matches[3].str());
            pos2 = stoi(matches[4].str());
          }
        }
      }
      else
        sequence += line;
    }
    if(!done)
      ig_file << "No probes found for gap " << lineCounter << std::endl;
  }
  fasta.close();
  pl_file.close();
  ip_file.close();
  ig_file.close();
}

int main(int argc, char **argv)
{
  int c;
  std::string genome = "hg19";
  std::string location = "chr8:133000000-133100000";
  std::string output_directory = "./output";
  while ((c = getopt (argc, argv, "g:l:o:")) != -1)
  {
    switch (c)
    {
      case 'g':
        genome = optarg;
        break;
      case 'l':
        location = optarg;
        break;
      case 'o':
        output_directory = optarg;
        break;
      case '?':
        return 1;
    }
  }

  const std::string restriction_enzyme = "DpnII";
  const std::string input_directory = ".";
  const std::string re_filename = input_directory + '/' + genome + '_' + restriction_enzyme + ".txt";
  const std::string tb_bin_filename = input_directory + "/twoBitToFa";
  const std::string tb_filename = input_directory + "/hg19.2bit";
  const std::string rs_filename = output_directory + "/restriction_sites.txt";
  const std::string tcb_filename = output_directory + "/tmp_candidates.bed";
  const std::string tcf_filename = output_directory + "/tmp_candidates.fa";
  const std::string pnp_filename = output_directory + "/probes_no_primers.txt";
  const std::string ip_filename = output_directory + "/issue_probes.txt";
  const std::string is_filename = output_directory + "/issue_restriction_sites.txt";
  const std::string pnps_filename = output_directory + "/probes_no_primers_sorted.txt";
  const std::string pnpno_filename = output_directory + "/probes_no_primers_no_overlaps.txt";
  const std::string gaps_filename = output_directory + "/gaps.txt";
  const std::string pnpnop2_filename = output_directory + "/probes_no_primers_no_overlaps_pass2.txt";
  const std::string ipp2_filename = output_directory + "/issue_probes_pass2.txt";
  const std::string ig_filename = output_directory + "/issue_gaps.txt";
  const std::string pnpnop2s_filename = output_directory + "/probes_no_primers_no_overlaps_pass2_sorted.txt";
  const std::string gapsp2_filename = output_directory + "/gaps_pass2.txt";
  const std::string pnpnop3_filename = output_directory + "/probes_no_primers_no_overlaps_pass3.txt";
  const std::string ipp3_filename = output_directory + "/issue_probes_pass3.txt";
  const std::string igp2_filename = output_directory + "/issue_gaps_pass2.txt";
  const std::string pnpnop3s_filename = output_directory + "/probes_no_primers_no_overlaps_pass3_sorted.txt";
  const std::string pwpno_filename = output_directory + "/probes_w_primers_no_overlaps.txt";
  const int max_length_from_rs = 80;
  const int max_length_from_rs2 = 110;
  const std::string forward_primer = "ATCGCACCAGCGTGT";
  const std::string reverse_primer = "CACTGCGGCTCCTCA";

  std::string chromosome;
  std::string chromosome_number;
  int start;
  int end;
  std::smatch matches;
  if(std::regex_search(location, matches, FASTA_REGEX))
  {
    chromosome = matches[1].str();
    chromosome_number = matches[2].str();
    start = stoi(matches[3].str());
    end = stoi(matches[4].str());
  }

  std::vector<int> restriction_sites;
  std::ifstream re_file(re_filename);
  std::ofstream rs_file(rs_filename);
  if(re_file.good() && rs_file.good())
  {
    std::string line;
    while(getline(re_file, line))
    {
      if(line.substr(0, chromosome_number.length()) == chromosome_number)
      {
        std::istringstream iss(line.substr(chromosome_number.length()));
        std::string field;
        while(iss >> field)
        {
          int restriction_site = stoi(field);
          if(restriction_site >= start && restriction_site <= end)
          {
            restriction_sites.push_back(restriction_site);
            rs_file << chromosome_number << '\t' << restriction_site << std::endl;
          }
          if(restriction_site > end)
            break;
        }
      }
    }
  }
  re_file.close();
  rs_file.close();

  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> probes;
  for(std::vector<int>::iterator it = restriction_sites.begin(); it != restriction_sites.end(); ++it)
  {
    std::ofstream tcb_file(tcb_filename);
    if(tcb_file.good())
      for(int i = 0; i <= max_length_from_rs; i++)
        tcb_file << chromosome << '\t' << (*it)-i-120 << '\t' << (*it)-i << '\t' << 1 << std::endl;
    tcb_file.close();
    system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
    ScanProbesPass1(tcf_filename, 10, 48, 84, 60, 72, pnp_filename, ip_filename, is_filename, "upstream", probes);
    tcb_file.open(tcb_filename);
    if(tcb_file.good())
      for(int i = 0; i <= max_length_from_rs; i++)
        tcb_file << chromosome << '\t' << (*it)+i << '\t' << (*it)+i+120 << '\t' << 1 << std::endl;
    tcb_file.close();
    system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
    ScanProbesPass1(tcf_filename, 10, 48, 84, 60, 72, pnp_filename, ip_filename, is_filename, "downstream", probes);
  }

  std::sort(probes.begin(), probes.end(), compareChromThenStart);
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator last = std::unique(probes.begin(), probes.end());
  probes.erase(last, probes.end());
  std::ofstream pnps_file(pnps_filename);
  if(pnps_file.good())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
      pnps_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
  }
  pnps_file.close();

  chromosome.clear();
  start = -1;
  end = -1;
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> probes_no_overlap;
  std::vector<std::tuple<std::string, int, int>> gaps;
  std::ofstream pnpno_file(pnpno_filename);
  std::ofstream gaps_file(gaps_filename, std::ofstream::app);
  if(pnpno_file.good() && gaps_file.good())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
    {
      if(chromosome.empty() || std::get<0>(*it) != chromosome)
      {
        probes_no_overlap.push_back(*it);
        pnpno_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
          ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
        chromosome = std::get<0>(*it);
        start = std::get<1>(*it);
        end = std::get<2>(*it);
      }
      else if(std::get<1>(*it) >= end)
      {
        probes_no_overlap.push_back(*it);
        pnpno_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
          ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
        if(std::get<1>(*it) - end + 10 > 120)
        {
          gaps.push_back(std::make_tuple(std::get<0>(*it), end - 5, std::get<1>(*it) + 5));
          gaps_file << std::get<0>(*it) << '\t' << end - 5 << '\t' << std::get<1>(*it) + 5 << '\t' << 1 << std::endl;
        }
        chromosome = std::get<0>(*it);
        start = std::get<1>(*it);
        end = std::get<2>(*it);
      }
    }
  }
  pnpno_file.close();
  gaps_file.close();
  probes = probes_no_overlap;

  system(("cp " + pnpno_filename + ' ' + pnpnop2_filename).c_str());

  std::unordered_set<int> closeRS;
  for(std::vector<int>::iterator it = restriction_sites.begin(); it != restriction_sites.end(); ++it)
      for (int i = -max_length_from_rs2; i<=max_length_from_rs2; i++)
        closeRS.insert((*it) + i);

  int l = 0;
  for(std::vector<std::tuple<std::string, int, int>>::iterator it = gaps.begin(); it != gaps.end(); ++it)
  {
    l++;
    std::string chromosome = std::get<0>(*it);
    int start = std::get<1>(*it);
    int end = std::get<2>(*it);
    std::ofstream tcb_file(tcb_filename);
    if(tcb_file.good())
      for(int i = start; i <= end-120; i++)
        if(closeRS.find(i) != closeRS.end() || closeRS.find(i + 120) != closeRS.end())
          tcb_file << chromosome << '\t' << i << '\t' << i + 120 << '\t' << 1 << std::endl;
    tcb_file.close();
    system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
    ScanProbesPass2_3(tcf_filename, 20, 48, 84, pnpnop2_filename, ipp2_filename, ig_filename, l, probes);
  }

  std::sort(probes.begin(), probes.end(), compareChromThenStart);
  last = std::unique(probes.begin(), probes.end());
  probes.erase(last, probes.end());
  std::ofstream pnpnop2s_file(pnpnop2s_filename);
  if(pnpnop2s_file.good())
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
      pnpnop2s_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
  pnpnop2s_file.close();

  std::map<std::string, int> x;
  for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
  {
    std::string key = std::get<0>(*it) + ' ' + std::to_string((int)((std::get<1>(*it) + 60) / 5000));
    std::map<std::string, int>::iterator match = x.find(key);
    if(match == x.end())
      x.insert(std::pair<std::string, int>(key, 1));
    else
      match->second++;
  }

  chromosome.clear();
  int pos1;
  int pos2;
  gaps.clear();
  std::ofstream gapsp2_file(gapsp2_filename);
  if(gapsp2_file.good())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
    {
      if(chromosome.empty())
      {
        chromosome = std::get<0>(*it);
        pos1 = std::get<1>(*it);
        pos2 = std::get<2>(*it);
      }
      else
      {
        int start = std::get<1>(*it);
        if (std::get<0>(*it)==chromosome && start-pos2+10>120 &&
          (x[chromosome + ' ' + std::to_string((int)(pos2/5000))]<5 ||
            x[chromosome + ' ' + std::to_string((int)(start/5000))]<5))
        {
          gaps.push_back(std::make_tuple(chromosome, pos2, start));
          gapsp2_file << chromosome << ' ' << pos2 << ' ' << start << std::endl;
        }
        chromosome = std::get<0>(*it);
        pos1 = std::get<1>(*it);
        pos2 = std::get<2>(*it);
      }
    }
  }
  gapsp2_file.close();

  system(("cp " + pnpnop2s_filename + ' ' + pnpnop3_filename).c_str());

  l = 0;
  for(std::vector<std::tuple<std::string, int, int>>::iterator it = gaps.begin(); it != gaps.end(); ++it)
  {
    l++;
    std::string chromosome = std::get<0>(*it);
    int start = std::get<1>(*it);
    int end = std::get<2>(*it);
    std::ofstream tcb_file(tcb_filename);
    if(tcb_file.good())
      for(int i = start; i <= end-120; i++)
        if(closeRS.find(i) != closeRS.end() || closeRS.find(i + 120) != closeRS.end())
          tcb_file << chromosome << '\t' << i << '\t' << i + 120 << '\t' << 1 << std::endl;
    tcb_file.close();
    system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
    ScanProbesPass2_3(tcf_filename, 25, 30, 96, pnpnop3_filename, ipp3_filename, igp2_filename, l, probes);
  }

  std::sort(probes.begin(), probes.end(), compareChromThenStart);
  last = std::unique(probes.begin(), probes.end());
  probes.erase(last, probes.end());
  std::ofstream pnpnop3s_file(pnpnop3s_filename);
  if(pnpnop3s_file.good())
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
      pnpnop3s_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
  pnpnop3s_file.close();

  std::ofstream pwpno_file(pwpno_filename);
  if(pwpno_file.good())
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
      pwpno_file << std::get<0>(*it) << '\t' << std::get<1>(*it) << '\t' << std::get<2>(*it) << '\t' <<
        forward_primer + std::get<3>(*it) + reverse_primer << '\t' << std::get<4>(*it) + forward_primer.length() + reverse_primer.length() << std::endl;
  pwpno_file.close();

  return 0;
}
