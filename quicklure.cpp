#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
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

void countRepeatsAndGC(const std::string sequence, int &repeats, int &GC)
{
  repeats = 0;
  int C = 0;
  int G = 0;
  for(int i = 0; i < sequence.length(); i++)
  {
    if(islower(sequence[i]))
      repeats++;
    if(sequence[i] == 'C')
      C++;
    if(sequence[i] == 'G')
      G++;
  }
  GC = G + C;
}

void ScanProbesPass1(std::vector<std::tuple<std::string, int, int>> probeCandidates, std::string roiSequence, int roiSequenceStart,
  int maximumRepeats, int relaxedMinimumGC, int relaxedMaximumGC, int stringentMinimumGC, int stringentMaximumGC,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueSitesFilename, std::string side,
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> &probes)
{
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream is_file(issueSitesFilename, std::ofstream::app);
  if(pl_file.good() && ip_file.good() && is_file.good())
  {
    std::tuple<std::string, int, int, std::string, int, int, int> backup_probe;
    bool done = false;
    for(std::vector<std::tuple<std::string, int, int>>::iterator it = probeCandidates.begin(); it != probeCandidates.end(); ++it)
    {
      std::string chromosome;
      int start;
      int end;
      std::tie(chromosome, start, end) = *it;
      std::string description = '>' + chromosome + ':' + std::to_string(start) + '-' + std::to_string(end);
      std::string sequence = roiSequence.substr(start - roiSequenceStart, end - start);
      int repeats;
      int GC;
      countRepeatsAndGC(sequence, repeats, GC);
      if(repeats <= maximumRepeats && GC >= relaxedMinimumGC && GC <= relaxedMaximumGC)
      {
        if(GC >= stringentMinimumGC && GC <= stringentMaximumGC)
        {
          probes.push_back(std::make_tuple(chromosome, start, end, sequence, sequence.length(), repeats, GC));
          pl_file << description << ' ' << sequence << ' ' << sequence.length() << ' ' << repeats << ' ' << GC << std::endl;
          done = true;
          break;
        }
        else if(std::get<0>(backup_probe).empty())
          backup_probe = std::make_tuple(chromosome, start, end, sequence, sequence.length(), repeats, GC);
      }
      else
        ip_file << description << ' ' << repeats << ' ' << GC << std::endl;
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
        std::vector<std::tuple<std::string, int, int>>::iterator first = probeCandidates.begin();
        std::string chromosome;
        int start;
        int end;
        std::tie(chromosome, start, end) = *first;
        is_file << "No " << side << " probes found for:" << chromosome << ' ' << (side == "upstream" ? end : start) << std::endl;
      }
    }
    pl_file.close();
    ip_file.close();
    is_file.close();
  }
}

void ScanProbesPass2_3(std::vector<std::tuple<std::string, int, int>> probeCandidates, std::string roiSequence, int roiSequenceStart,
  int maximumRepeats, int minimumGC, int maximumGC,
  std::string probeListFilename, std::string issueProbesFilename, std::string issueGapsFilename, int gapNumber,
  std::vector<std::tuple<std::string, int, int, std::string, int, int, int>> &probes)
{
  std::ofstream pl_file(probeListFilename, std::ofstream::app);
  std::ofstream ip_file(issueProbesFilename, std::ofstream::app);
  std::ofstream ig_file(issueGapsFilename, std::ofstream::app);
  if(pl_file.good() && ip_file.good() && ig_file.good())
  {
    std::string chromosome;
    int pos1;
    int pos2;
    bool done = false;
    int end = 0;
    bool update = true;
    for(std::vector<std::tuple<std::string, int, int>>::iterator it = probeCandidates.begin(); it != probeCandidates.end(); ++it)
    {
      if(update)
      {
        std::tie(chromosome, pos1, pos2) = *it;
        update = false;
      }
      if(pos1 > end)
      {
        std::string sequence = roiSequence.substr(pos1 - roiSequenceStart, pos2 - pos1);
        int repeats;
        int GC;
        countRepeatsAndGC(sequence, repeats, GC);
        if(repeats <= maximumRepeats && GC >= minimumGC && GC <= maximumGC)
        {
          probes.push_back(std::make_tuple(chromosome, pos1, pos2, sequence, sequence.length(), repeats, GC));
          pl_file << chromosome << ' ' << pos1 << ' ' << pos2 << ' ' << sequence << ' ' << sequence.length() << ' ' << repeats << ' ' << GC << std::endl;
          end = pos2;
          done = true;
        }
        else
          ip_file << '>' << chromosome << ':' << pos1 << '-' << pos2 << ' ' << repeats << ' ' << GC << std::endl;
        update = true;
      }
    }
    if(!done)
      ig_file << "No probes found for gap " << gapNumber << std::endl;
  }
  pl_file.close();
  ip_file.close();
  ig_file.close();
}

int main(int argc, char **argv)
{
  int c;
  std::string genome = "hg19";
  std::string restriction_enzyme = "DpnII";
  std::string location = "chr8:133000000-133100000";
  std::string output_directory = "./output";
  while((c = getopt(argc, argv, "g:e:l:o:h")) != -1)
  {
    switch (c)
    {
      case 'g':
        genome = optarg;
        break;
      case 'e':
        restriction_enzyme = optarg;
        break;
      case 'l':
        location = optarg;
        break;
      case 'o':
        output_directory = optarg;
        break;
      case 'h':
        std::cout << "Usage: " << argv[0] <<
          " [-g genome ID (default: " << genome << ")]" <<
          " [-e restriction enzyme (default: " << restriction_enzyme << ")]" <<
          " [-l location (default: " << location << ")]" <<
          " [-o output directory (default: " << output_directory << ") ]" <<
          " [-h]" << std::endl;
        return 0;
    }
  }

  const std::string input_directory = ".";
  const std::string re_filename = input_directory + '/' + genome + '_' + restriction_enzyme + ".txt";
  const std::string tb_bin_filename = input_directory + "/twoBitToFa";
  const std::string tb_filename = input_directory + '/' + genome + ".2bit";
  const std::string rs_filename = output_directory + "/restriction_sites.txt";
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
  const int probe_length = 120;
  const int max_length_from_rs = 80;
  const int max_length_from_rs2 = 110;
  const int pass1_maximum_repeats = 10;
  const int pass1_relaxed_minimum_GC = 48;
  const int pass1_relaxed_maximum_GC = 84;
  const int pass1_stringent_minimum_GC = 60;
  const int pass1_stringent_maximum_GC = 72;
  const int pass2_maximum_repeats = 20;
  const int pass2_minimum_GC = 48;
  const int pass2_maximum_GC = 84;
  const int pass3_maximum_repeats = 25;
  const int pass3_minimum_GC = 30;
  const int pass3_maximum_GC = 96;
  const std::string forward_primer = "ATCGCACCAGCGTGT";
  const std::string reverse_primer = "CACTGCGGCTCCTCA";

  std::string chromosome;
  std::string chromosome_number;
  int roi_start;
  int roi_end;
  std::smatch matches;
  if(std::regex_search(location, matches, FASTA_REGEX))
  {
    chromosome = matches[1].str();
    chromosome_number = matches[2].str();
    roi_start = stoi(matches[3].str());
    roi_end = stoi(matches[4].str());
  }

  const std::string tb_extension = ".2bit";
  int index = tb_filename.rfind(tb_extension);
  std::string tf_filename = tb_filename;
  tf_filename.replace(index, tb_extension.length(), ".fa");
  system((tb_bin_filename + ' ' + tb_filename + ':' + chromosome + ':' + std::to_string(roi_start) + '-' + std::to_string(roi_end) + ' ' + tf_filename).c_str());
  std::string roi_sequence;
  std::ifstream tf_file(tf_filename);
  if(tf_file.good())
  {
    std::string line;
    while(getline(tf_file, line))
      if(line[0] != '>')
        roi_sequence += line;
  }
  tf_file.close();

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
          if(restriction_site >= roi_start && restriction_site <= roi_end)
          {
            restriction_sites.push_back(restriction_site);
            rs_file << chromosome_number << '\t' << restriction_site << std::endl;
          }
          if(restriction_site > roi_end)
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
    std::vector<std::tuple<std::string, int, int>> probe_candidates;
    for(int i = 0; i <= max_length_from_rs; i++)
      probe_candidates.push_back(std::make_tuple(chromosome, (*it) - i - probe_length, (*it) - i));
    ScanProbesPass1(probe_candidates, roi_sequence, roi_start,
      pass1_maximum_repeats, pass1_relaxed_minimum_GC, pass1_relaxed_maximum_GC, pass1_stringent_minimum_GC, pass1_stringent_maximum_GC,
      pnp_filename, ip_filename, is_filename, "upstream", probes);
    probe_candidates.clear();
    for(int i = 0; i <= max_length_from_rs; i++)
      probe_candidates.push_back(std::make_tuple(chromosome, (*it) + i, (*it) + i + probe_length));
    ScanProbesPass1(probe_candidates, roi_sequence, roi_start,
      pass1_maximum_repeats, pass1_relaxed_minimum_GC, pass1_relaxed_maximum_GC, pass1_stringent_minimum_GC, pass1_stringent_maximum_GC,
      pnp_filename, ip_filename, is_filename, "downstream", probes);
  }

  std::sort(probes.begin(), probes.end());
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
  int start;
  int end;
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
        std::tie(chromosome, start, end, std::ignore, std::ignore, std::ignore, std::ignore) = *it;
      }
      else if(std::get<1>(*it) >= end)
      {
        probes_no_overlap.push_back(*it);
        pnpno_file << std::get<0>(*it) << ' ' << std::get<1>(*it) << ' ' << std::get<2>(*it) << ' ' << std::get<3>(*it) <<
          ' ' << std::get<4>(*it) << ' ' << std::get<5>(*it) << ' ' << std::get<6>(*it) << std::endl;
        if(std::get<1>(*it) - end + 10 > probe_length)
        {
          gaps.push_back(std::make_tuple(std::get<0>(*it), end - 5, std::get<1>(*it) + 5));
          gaps_file << std::get<0>(*it) << '\t' << end - 5 << '\t' << std::get<1>(*it) + 5 << '\t' << 1 << std::endl;
        }
        std::tie(chromosome, start, end, std::ignore, std::ignore, std::ignore, std::ignore) = *it;
      }
    }
  }
  pnpno_file.close();
  gaps_file.close();
  probes = probes_no_overlap;

  system(("cp " + pnpno_filename + ' ' + pnpnop2_filename).c_str());

  std::unordered_set<int> closeRS;
  for(std::vector<int>::iterator it = restriction_sites.begin(); it != restriction_sites.end(); ++it)
    for(int i = -max_length_from_rs2; i <= max_length_from_rs2; i++)
      closeRS.insert((*it) + i);

  int l = 0;
  for(std::vector<std::tuple<std::string, int, int>>::iterator it = gaps.begin(); it != gaps.end(); ++it)
  {
    l++;
    std::string chromosome;
    int start;
    int end;
    std::tie(chromosome, start, end) = *it;
    std::vector<std::tuple<std::string, int, int>> probe_candidates;
    for(int i = start; i <= end - probe_length; i++)
      if(closeRS.find(i) != closeRS.end() || closeRS.find(i + probe_length) != closeRS.end())
        probe_candidates.push_back(std::make_tuple(chromosome, i, i + probe_length));
    ScanProbesPass2_3(probe_candidates, roi_sequence, roi_start,
      pass2_maximum_repeats, pass2_minimum_GC, pass2_maximum_GC,
      pnpnop2_filename, ipp2_filename, ig_filename, l, probes);
  }

  std::sort(probes.begin(), probes.end());
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
        std::tie(chromosome, pos1, pos2, std::ignore, std::ignore, std::ignore, std::ignore) = *it;
      else
      {
        int start = std::get<1>(*it);
        if (std::get<0>(*it) == chromosome && start - pos2 + 10 > probe_length &&
          (x[chromosome + ' ' + std::to_string((int)(pos2 / 5000))] < 5 ||
          x[chromosome + ' ' + std::to_string((int)(start / 5000))] < 5))
        {
          gaps.push_back(std::make_tuple(chromosome, pos2, start));
          gapsp2_file << chromosome << ' ' << pos2 << ' ' << start << std::endl;
        }
        std::tie(chromosome, pos1, pos2, std::ignore, std::ignore, std::ignore, std::ignore) = *it;
      }
    }
  }
  gapsp2_file.close();

  system(("cp " + pnpnop2s_filename + ' ' + pnpnop3_filename).c_str());

  l = 0;
  for(std::vector<std::tuple<std::string, int, int>>::iterator it = gaps.begin(); it != gaps.end(); ++it)
  {
    l++;
    std::string chromosome;
    int start;
    int end;
    std::tie(chromosome, start, end) = *it;
    std::vector<std::tuple<std::string, int, int>> probe_candidates;
    for(int i = start; i <= end - probe_length; i++)
      if(closeRS.find(i) != closeRS.end() || closeRS.find(i + probe_length) != closeRS.end())
        probe_candidates.push_back(std::make_tuple(chromosome, i, i + probe_length));
    ScanProbesPass2_3(probe_candidates, roi_sequence, roi_start,
      pass3_maximum_repeats, pass3_minimum_GC, pass3_maximum_GC,
      pnpnop3_filename, ipp3_filename, igp2_filename, l, probes);
  }

  std::sort(probes.begin(), probes.end());
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
