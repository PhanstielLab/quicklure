#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <regex>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <vector>

void ScanProbesPass1(std::string candidatesFastaFilename, int repeatThresh, int GCThresh1, int GCThresh2, int GCThresh3, int GCThresh4,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueSitesFilename, std::string side)
{
  std::ifstream fasta(candidatesFastaFilename);
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream is_file(issueSitesFilename, std::ofstream::app);
  if(fasta.is_open() && pl_file.is_open() && ip_file.is_open() && is_file.is_open())
  {
    std::string line;
    std::string description;
    int entries = 0;
    std::string sequence;
    std::string first_description;
    std::string backup;
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
            if(c1_file.is_open())
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
            if(GC>=GCThresh3&&GC<=GCThresh4)
            {
              std::ofstream s("something.txt", std::ofstream::app);
              if(s.is_open())
                s << sequence << std::endl;
              s.close();
              pl_file << description << ' ' << sequence << ' ' << sequence.length() << ' ' << counter << ' ' << GC << std::endl;
              done = true;
              break;
            }
            else if(backup.empty())
              backup = description + ' ' + sequence + ' ' + std::to_string(sequence.length()) + ' ' + std::to_string(counter) + ' ' + std::to_string(GC);
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
      if(!backup.empty())
        pl_file << backup << std::endl;
      else
      {
        std::regex rgx(">(chr.+):(\\d+)-(\\d+)");
        std::smatch matches;
        std::regex_search(first_description, matches, rgx);
        std::string chromosome = matches[1].str();
        int start = stoi(matches[2].str());
        int end = stoi(matches[3].str());
        is_file << "No " << side << " probes found for:" << chromosome << ' ' << (side == "upstream" ? end : start) << std::endl;
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

bool compareProbes(const std::tuple<std::string, int, int, std::string, int, int, int> &left,
  const std::tuple<std::string, int, int, std::string, int, int, int> &right)
{
  return std::get<0>(left) == std::get<0>(right) && std::get<1>(left) == std::get<1>(right) &&
    std::get<2>(left) == std::get<2>(right) && std::get<3>(left) == std::get<3>(right) &&
    std::get<4>(left) == std::get<4>(right) && std::get<5>(left) == std::get<5>(right) &&
    std::get<6>(left) == std::get<6>(right);
}

void ScanProbesPass2_3(std::string candidatesFastaFilename, int repeatThresh, int GCThresh1, int GCThresh2,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueGapsFilename, int lineCounter)
{
  std::ifstream fasta(candidatesFastaFilename);
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream ig_file(issueGapsFilename, std::ofstream::app);
  if(fasta.is_open() && pl_file.is_open() && ip_file.is_open() && ig_file.is_open())
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
              pl_file << chromosome << ' ' << pos1 << ' ' << pos2 << ' ' << sequence << ' ' << sequence.length() << ' ' << counter << ' ' << GC << std::endl;
              end = pos2;
              done = true;
            }
            else
              ip_file << description << ' ' << counter << ' ' << GC << std::endl;
            std::regex rgx(">(.+):(\\d+)-(\\d+)");
            std::smatch matches;
            std::regex_search(description, matches, rgx);
            chromosome = matches[1].str();
            pos1 = stoi(matches[2].str());
            pos2 = stoi(matches[3].str());
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
          std::regex rgx(">(.+):(\\d+)-(\\d+)");
          std::smatch matches;
          std::regex_search(description, matches, rgx);
          chromosome = matches[1].str();
          pos1 = stoi(matches[2].str());
          pos2 = stoi(matches[3].str());
        }
      }
      else
        sequence += line;
    }
    if(!done)
      ig_file << "No probes found for gap " << lineCounter << std::endl;
    fasta.close();
    pl_file.close();
    ip_file.close();
    ig_file.close();
  }
}

int main()
{
  std::string genome = "hg19";
  std::string restriction_enzyme = "DpnII";
  std::string re_filename = "/Users/cwenger/Desktop/Doug/lure++/" + genome + '_' + restriction_enzyme + ".txt";
  std::string rs_filename = "/Users/cwenger/Desktop/Doug/lure++/output/restriction_sites.txt";
  std::string tcb_filename = "/Users/cwenger/Desktop/Doug/lure++/output/tmp_candidates.bed";
  std::string tcf_filename = "/Users/cwenger/Desktop/Doug/lure++/output/tmp_candidates.fa";
  std::string tb_bin_filename = "/Users/cwenger/Desktop/Doug/lure++/eric_scripts/twoBitToFa";
  std::string tb_filename = "/Users/cwenger/Desktop/Doug/lure++/eric_scripts/hg19.2bit";
  std::string pnp_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers.txt";
  std::string ip_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_probes.txt";
  std::string is_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_restriction_sites.txt";
  std::string pnps_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_sorted.txt";
  std::string pnpno_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_no_overlaps.txt";
  std::string gaps_filename = "/Users/cwenger/Desktop/Doug/lure++/output/gaps.txt";
  std::string pnpnop2_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_no_overlaps_pass2.txt";
  std::string ipp2_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_probes_pass2.txt";
  std::string ig_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_gaps.txt";
  std::string pnpnop2s_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_no_overlaps_pass2_sorted.txt";
  std::string gapsp2_filename = "/Users/cwenger/Desktop/Doug/lure++/output/gaps_pass2.txt";
  std::string pnpnop3_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_no_overlaps_pass3.txt";
  std::string ipp3_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_probes_pass3.txt";
  std::string igp2_filename = "/Users/cwenger/Desktop/Doug/lure++/output/issue_gaps_pass2.txt";
  std::string pnpnop3s_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_no_primers_no_overlaps_pass3_sorted.txt";
  std::string pwpno_filename = "/Users/cwenger/Desktop/Doug/lure++/output/probes_w_primers_no_overlaps.txt";
  int max_length_from_rs = 80;
  int max_length_from_rs2 = 110;
  std::string forward_primer = "ATCGCACCAGCGTGT";
  std::string reverse_primer = "CACTGCGGCTCCTCA";

  std::string location = "chr8:133000000-134000000";
  std::regex rgx("(chr.+):(\\d+)-(\\d+)");
  std::smatch matches;
  std::regex_search(location, matches, rgx);
  std::string chromosome = matches[1].str();
  std::string chromosome_number = chromosome.substr(3);
  int start = stoi(matches[2].str());
  int end = stoi(matches[3].str());

  std::ifstream re_file(re_filename);
  std::ofstream rs_file_out(rs_filename);
  if(re_file.is_open() && rs_file_out.is_open())
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
          int residue = stoi(field);
          if(residue >= start && residue <= end)
            rs_file_out << chromosome_number << '\t' << residue << std::endl;
          if(residue > end)
            break;
        }
      }
    }
    re_file.close();
    rs_file_out.close();
  }

  std::ifstream rs_file_in(rs_filename);
  if(rs_file_in.is_open())
  {
    std::string line;
    while(getline(rs_file_in, line))
    {
      int tab = line.find('\t');
      int res = stoi(line.substr(tab + 1));
      std::ofstream tcb_file(tcb_filename);
      if(tcb_file.is_open())
      {
        for(int i = 0; i <= max_length_from_rs; i++)
          tcb_file << chromosome << '\t' << res-i-120 << '\t' << res-i << '\t' << 1 << std::endl;
        tcb_file.close();
        system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
        ScanProbesPass1(tcf_filename, 10, 48, 84, 60, 72, pnp_filename, ip_filename, is_filename, "upstream");
      }
      tcb_file.open(tcb_filename);
      if(tcb_file.is_open())
      {
        for(int i = 0; i <= max_length_from_rs; i++)
          tcb_file << chromosome << '\t' << res+i << '\t' << res+i+120 << '\t' << 1 << std::endl;
        tcb_file.close();
        system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
        ScanProbesPass1(tcf_filename, 10, 48, 84, 60, 72, pnp_filename, ip_filename, is_filename, "downstream");
      }
    }

    rs_file_in.close();
  }

  std::ifstream pnp_file(pnp_filename);
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> probes;
  if(pnp_file.is_open())
  {
    std::string line;
    while(getline(pnp_file, line))
    {
      std::regex rgx(">(.+):(\\d+)-(\\d+) (.+) (\\d+) (\\d+) (\\d+)");
      std::smatch matches;
      std::regex_search(line, matches, rgx);
      std::string chromosome = matches[1].str();
      int start = stoi(matches[2].str());
      int end = stoi(matches[3].str());
      std::string sequence = matches[4].str();
      int length = stoi(matches[5].str());
      int counter = stoi(matches[6].str());
      int GC = stoi(matches[7].str());
      probes.push_back(std::tuple<std::string, int, int, std::string, int, int, int>(chromosome, start, end, sequence, length, counter, GC));
    }
    pnp_file.close();
  }
  std::sort(probes.begin(), probes.end(), compareChromThenStart);
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it2 = std::unique(probes.begin(), probes.end(), compareProbes);
  probes.resize(std::distance(probes.begin(), it2));
  std::ofstream pnps_file(pnps_filename);
  if(pnps_file.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
      pnps_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
    pnps_file.close();
  }

  chromosome.clear();
  start = -1;
  end = -1;
  std::ofstream pnpno_file(pnpno_filename);
  std::ofstream gaps_file_out(gaps_filename, std::ofstream::app);
  if(pnpno_file.is_open() && gaps_file_out.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes.begin(); it != probes.end(); ++it)
    {
      if(chromosome.empty() || std::get<0>(*it) != chromosome)
      {
        pnpno_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
          ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
        chromosome = std::get<0>(*it);
        start = std::get<1>(*it);
        end = std::get<2>(*it);
      }
      else if(std::get<1>(*it) >= end)
      {
        pnpno_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
          ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
        if(std::get<1>(*it) - end + 10 > 120)
          gaps_file_out << std::get<0>(*it) << '\t' << end - 5 << '\t' << std::get<1>(*it) + 5 << '\t' << 1 << std::endl;
        chromosome = std::get<0>(*it);
        start = std::get<1>(*it);
        end = std::get<2>(*it);
      }
    }
    pnpno_file.close();
    gaps_file_out.close();
  }

  system(("cp " + pnpno_filename + ' ' + pnpnop2_filename).c_str());

  std::unordered_set<int> closeRS;
  rs_file_in.open(rs_filename);
  if(rs_file_in.is_open())
  {
    std::string line;
    while(getline(rs_file_in, line))
    {
      int tab = line.find('\t');
      int res = stoi(line.substr(tab + 1));
      for (int i = -max_length_from_rs2; i<=max_length_from_rs2; i++)
        closeRS.insert(res + i);
    }

    rs_file_in.close();
  }

  std::ifstream gaps_file_in(gaps_filename);
  int l = 0;
  if(gaps_file_in.is_open())
  {
    std::string line;
    while(getline(gaps_file_in, line))
    {
      l++;
      std::regex rgx("(.+)\\t(\\d+)\\t(\\d+)\\t(\\d+)");
      std::smatch matches;
      std::regex_search(line, matches, rgx);
      std::string chromosome = matches[1].str();
      int start = stoi(matches[2].str());
      int end = stoi(matches[3].str());
      std::ofstream tcb_file(tcb_filename);
      if(tcb_file.is_open())
      {
        for(int i = start; i <= end-120; i++)
          if(closeRS.find(i) != closeRS.end() || closeRS.find(i + 120) != closeRS.end())
            tcb_file << chromosome << '\t' << i << '\t' << i + 120 << '\t' << 1 << std::endl;
        tcb_file.close();
      }
      system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
      ScanProbesPass2_3(tcf_filename, 20, 48, 84, pnpnop2_filename, ipp2_filename, ig_filename, l);
    }

    gaps_file_in.close();
  }

  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> probes2;
  std::ifstream pnpnop2_file(pnpnop2_filename);
  if(pnpnop2_file.is_open())
  {
    std::string line;
    while(getline(pnpnop2_file, line))
    {
      std::regex rgx("(.+) (\\d+) (\\d+) (.+) (\\d+) (\\d+) (\\d+)");
      std::smatch matches;
      std::regex_search(line, matches, rgx);
      std::string chromosome = matches[1].str();
      int start = stoi(matches[2].str());
      int end = stoi(matches[3].str());
      std::string sequence = matches[4].str();
      int length = stoi(matches[5].str());
      int counter = stoi(matches[6].str());
      int GC = stoi(matches[7].str());
      probes2.push_back(std::make_tuple(chromosome, start, end, sequence, length, counter, GC));
    }

    pnpnop2_file.close();
  }
  std::sort(probes2.begin(), probes2.end(), compareChromThenStart);
  it2 = std::unique(probes2.begin(), probes2.end(), compareProbes);
  probes2.resize(std::distance(probes2.begin(), it2));
  std::ofstream pnpnop2s_file(pnpnop2s_filename);
  if(pnpnop2s_file.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes2.begin(); it != probes2.end(); ++it)
      pnpnop2s_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
    pnpnop2s_file.close();
  }

  std::map<std::string, int> x;
  for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes2.begin(); it != probes2.end(); ++it)
  {
    std::string key = std::get<0>(*it) + ' ' + std::to_string((int)((std::get<1>(*it) + 60) / 5000));
    std::map<std::string, int>::iterator it2 = x.find(key);
    if(it2 == x.end())
      x.insert(std::pair<std::string, int>(key, 1));
    else
      it2->second++;
  }

  chromosome.clear();
  int pos1;
  int pos2;
  std::ofstream gapsp2_file_out(gapsp2_filename);
  if(gapsp2_file_out.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes2.begin(); it != probes2.end(); ++it)
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
          gapsp2_file_out << chromosome << ' ' << pos2 << ' ' << start << std::endl;
        chromosome = std::get<0>(*it);
        pos1 = std::get<1>(*it);
        pos2 = std::get<2>(*it);
      }
    }

    gapsp2_file_out.close();
  }

  system(("cp " + pnpnop2s_filename + ' ' + pnpnop3_filename).c_str());

  // can't we just use closeRS from earlier?
  closeRS.clear();
  rs_file_in.open(rs_filename);
  if(rs_file_in.is_open())
  {
    std::string line;
    while(getline(rs_file_in, line))
    {
      int tab = line.find('\t');
      int res = stoi(line.substr(tab + 1));
      for (int i = -max_length_from_rs2; i<=max_length_from_rs2; i++)
        closeRS.insert(res + i);
    }

    rs_file_in.close();
  }

  std::ifstream gapsp2_file_in(gapsp2_filename);
  l = 0;
  if(gapsp2_file_in.is_open())
  {
    std::string line;
    while(getline(gapsp2_file_in, line))
    {
      l++;
      std::regex rgx("(.+) (\\d+) (\\d+)");
      std::smatch matches;
      std::regex_search(line, matches, rgx);
      std::string chromosome = matches[1].str();
      int start = stoi(matches[2].str());
      int end = stoi(matches[3].str());
      std::ofstream tcb_file(tcb_filename);
      if(tcb_file.is_open())
      {
        for(int i = start; i <= end-120; i++)
          if(closeRS.find(i) != closeRS.end() || closeRS.find(i + 120) != closeRS.end())
            tcb_file << chromosome << '\t' << i << '\t' << i + 120 << '\t' << 1 << std::endl;
        tcb_file.close();
      }
      system((tb_bin_filename + ' ' + tb_filename + ' ' + tcf_filename + " -bed=" + tcb_filename + " -bedPos").c_str());
      ScanProbesPass2_3(tcf_filename, 25, 30, 96, pnpnop3_filename, ipp3_filename, igp2_filename, l);
    }

    gapsp2_file_in.close();
  }

  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> probes3;
  std::ifstream pnpnop3_file(pnpnop3_filename);
  if(pnpnop3_file.is_open())
  {
    std::string line;
    while(getline(pnpnop3_file, line))
    {
      std::regex rgx("(.+) (\\d+) (\\d+) (.+) (\\d+) (\\d+) (\\d+)");
      std::smatch matches;
      std::regex_search(line, matches, rgx);
      std::string chromosome = matches[1].str();
      int start = stoi(matches[2].str());
      int end = stoi(matches[3].str());
      std::string sequence = matches[4].str();
      int length = stoi(matches[5].str());
      int counter = stoi(matches[6].str());
      int GC = stoi(matches[7].str());
      probes3.push_back(std::make_tuple(chromosome, start, end, sequence, length, counter, GC));
    }

    pnpnop3_file.close();
  }
  std::sort(probes3.begin(), probes3.end(), compareChromThenStart);
  it2 = std::unique(probes3.begin(), probes3.end(), compareProbes);
  probes3.resize(std::distance(probes3.begin(), it2));
  std::ofstream pnpnop3s_file(pnpnop3s_filename);
  if(pnpnop3s_file.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes3.begin(); it != probes3.end(); ++it)
      pnpnop3s_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
        ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
    pnpnop3s_file.close();
  }

  std::ofstream pwpno_file(pwpno_filename);
  if(pwpno_file.is_open())
  {
    for(std::vector<std::tuple<std::string, int, int, std::string, int, int, int>>::iterator it = probes3.begin(); it != probes3.end(); ++it)
      pwpno_file << std::get<0>(*it) << '\t' << std::get<1>(*it) << '\t' << std::get<2>(*it) << '\t' <<
        forward_primer + std::get<3>(*it) + reverse_primer << '\t' << std::get<4>(*it) + forward_primer.length() + reverse_primer.length() << std::endl;
    pwpno_file.close();
  }

  return 0;
}
