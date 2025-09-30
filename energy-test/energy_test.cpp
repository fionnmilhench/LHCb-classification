#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <args.hxx>

double delta = 0.5;
double divisor = 2 / (256.0 * 2 * delta * delta);

#define EXACT 1

#if EXACT
double my_exp(double x) {
    return std::exp(-x/(delta*delta));
}
#else
double my_exp(double x) {
    x = std::max(1 - x * divisor, 0.0);
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return x;
}
#endif

struct Event {
  double s12;
  double s13;
  double s24;
  double s34;
  double s134;
  double half_mag_squared;
};

const std::vector<Event> read_file(const std::string filename, const size_t n_events) {
  std::fstream file(filename, std::ios_base::in);
  if (!file)
    throw std::runtime_error("Error opening file " + filename);

  std::vector<Event> events;
  events.reserve(std::min((size_t) 5000000, n_events));

  std::string line;
  while (std::getline(file, line) && events.size() < n_events) {
    std::istringstream iss(line);
    Event e;
    iss >> e.s12 >> e.s13 >> e.s24 >> e.s34 >> e.s134;
    if (iss.fail())
      throw std::runtime_error("Error reading line " + std::to_string(events.size()+1) + " in "  + filename);
    e.half_mag_squared = 0.5 * (e.s12*e.s12 + e.s13*e.s13 + e.s24*e.s24 + e.s34*e.s34 + e.s134*e.s134);
    events.push_back(e);
  }
  return events;
}

double compute_distance(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool upper_only) {
  double total = 0;
  #pragma omp parallel for reduction(+:total) schedule(static, 1)
  for (size_t i=0; i < data_1.size(); ++i) {
    auto event_1 = data_1[i];
    for (size_t j=(upper_only ? i+1 : 0); j < data_2.size(); ++j) {
      auto event_2 = data_2[j];
      double distance_squared = event_1.half_mag_squared + event_2.half_mag_squared - (
        event_1.s13 * event_2.s13 + event_1.s12 * event_2.s12 +
        event_1.s24 * event_2.s24 + event_1.s34 * event_2.s34 +
        event_1.s134 * event_2.s134);
      total += my_exp(distance_squared);
    }
  }
  return total;
}

double compute_statistic(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool debug = false) {
  double dist_11 = compute_distance(data_1, data_1, true);
  dist_11 /= data_1.size()*(data_1.size()-1);
  if (debug)
    std::cout << "    dist_11 = " << dist_11 << std::endl;

  double dist_22 = compute_distance(data_2, data_2, true);
  dist_22 /= data_2.size()*(data_2.size()-1);
  if (debug)
    std::cout << "    dist_22 = " << dist_22 << std::endl;

  double dist_12 = compute_distance(data_1, data_2, false);
  dist_12 /= data_1.size()*data_2.size();
  if (debug)
    std::cout << "    dist_12 = " << dist_12 << std::endl;

  return dist_11 + dist_22 - dist_12;
}

std::vector<double> compute_statistic_contributions(const std::vector<Event> &data_1, const std::vector<Event> &data_2, const bool debug = false) {
  // Create vector to contain individual test statistic contributions (Ti values)
  std::vector<double> statistic_contributions;
  statistic_contributions.reserve(data_1.size());
  
  // Calculate Ti values for all events in data_1
  for (size_t i=0; i < data_1.size(); ++i) {
    // Create length-1 vector to contain the current working event
    std::vector<Event> event;
    event.push_back(data_1[i]);
    
    double dist_e1 = compute_distance(event, data_1, false) - 1; // Subtract 1 to cancel the event's distance to itself
    dist_e1 /= 2 * data_1.size() * (data_1.size()-1);
    
    double dist_e2 = compute_distance(event, data_2, false);
    dist_e2 /= 2 * data_1.size() * data_2.size();
    
    statistic_contributions.push_back(dist_e1 - dist_e2);
  }
  
  // Ensure that vector length matches
  if (data_1.size() != statistic_contributions.size()) throw "Mismatch between number of Ti values and number of events";
  
  return statistic_contributions;
}

int run_energy_test(int argc, char *argv[]) {
  #if EXACT
    std::cout << "Running in exact mode" << std::endl;
  #else
    std::cout << "Running in approximate mode" << std::endl;
  #endif

  // Parse the command line arguments
  args::ArgumentParser parser("CPU based energy test");
  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  args::Flag calculate_ti(parser, "calculate ti", "Calculate individual contributions to test statistic", {"calculate-ti"});
  args::Flag permutations_only(parser, "permutations only", "Only calculate permutations", {"permutations-only"});
  args::Flag output_write(parser, "output write", "write output Tvalues", {"output-write"});
  args::ValueFlag<size_t> n_permutations(parser, "n_permutations", "Number of permutations to run", {"n-permutations"});
  args::ValueFlag<size_t> max_events_1(parser, "max events 1", "Maximum number of events to use from dataset 1", {"max-events-1"});
  args::ValueFlag<size_t> max_events_2(parser, "max events 2", "Maximum number of events to use from dataset 2", {"max-events-2"});
  args::ValueFlag<size_t> max_events(parser, "max events", "Max number of events in each dataset", {"max-events"});
  args::ValueFlag<size_t> seed(parser, "seed", "seed for permutations", {"seed"});
  args::ValueFlag<size_t> max_permutation_events_1(parser, "max permutation events 1", "Max number of events in dataset 1 for permutations",
                                                   {"max-permutation-events-1"});
  args::ValueFlag<double> delta_value(parser, "delta value", "delta_value", {"delta-value"});
  args::ValueFlag<std::string> ti_output_fn_1(parser, "ti output filename 1", "Output filename for individual contributions to test statistic from dataset 1", {"ti-output-fn-1"});
  args::ValueFlag<std::string> ti_output_fn_2(parser, "ti output filename 2", "Output filename for individual contributions to test statistic from dataset 2", {"ti-output-fn-2"});
  args::ValueFlag<std::string> permutation_ti_minmax_output_fn(parser, "permutation ti min-max filename", "Output filename for the minimum and maximum Ti values from permutations", {"permutation-ti-minmax-output-fn"});
  args::Positional<std::string> filename_1(parser, "dataset 1", "Filename for the first dataset");
  args::Positional<std::string> filename_2(parser, "dataset 2", "Filename for the second dataset");
  args::Positional<std::string> permutation_output_fn(parser, "permutation output filename", "Output filename for the permutation test statistics", {"permutation-output-fn"});


   
  try {
    parser.ParseCLI(argc, argv);
    if (!filename_1 || !filename_2)
      throw args::ParseError("Two dataset filenames must be given");
    if ((max_events_1 || max_events_2) && max_events)
      throw args::ParseError("--max-events cannot be used with --max-events-1 or --max-events-2");
    if (calculate_ti && max_permutation_events_1)
      throw args::ParseError("--calculate-ti cannot be used with --max-permutation-events-1");
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  //set delta
  if (delta_value){
    delta = args::get(delta_value);
  }
  std::cout<<"Distance parameter set to "<< delta <<std::endl;
  divisor = 2 / (256.0 * 2 * delta * delta);
 
    
  
  // Parse the maximum number of events to use
  size_t data_1_limit = std::numeric_limits<size_t>::max();
  size_t data_2_limit = std::numeric_limits<size_t>::max();

  if (max_events) {
    data_1_limit = args::get(max_events);
    data_2_limit = args::get(max_events);
  } else {
    if (max_events_1)
      data_1_limit = args::get(max_events_1);
    if (max_events_2)
      data_2_limit = args::get(max_events_2);
  }

  // Load the data
  const auto dataset_1 = read_file(args::get(filename_1), data_1_limit);
  std::cout << "Dataset 1 size is " << dataset_1.size() << std::endl;
  const auto dataset_2 = read_file(args::get(filename_2), data_2_limit);
  std::cout << "Dataset 2 size is " << dataset_2.size() << std::endl;
  std::cout << std::endl;
  
  double real_test_statistic = -999;

  if(!permutations_only){
    // Compute the individual test statistic contributions (Ti values)
    if (calculate_ti) {
      // Get Ti vectors
      std::cout << "Calculating contributions of individual events to test statistic..." << std::endl;
      const std::vector<double> tis_1 = compute_statistic_contributions(dataset_1, dataset_2);
      std::cout << "  Finished dataset 1" << std::endl;
      const std::vector<double> tis_2 = compute_statistic_contributions(dataset_2, dataset_1);
      std::cout << "  Finished dataset 2" << std::endl << std::endl;
      
      // Display full test statistic
      double total = 0;
      for (size_t i=0; i < tis_1.size(); ++i) total += tis_1[i];
      for (size_t i=0; i < tis_2.size(); ++i) total += tis_2[i];
      real_test_statistic = total;
      std::cout << "Test statistic for nominal dataset:" << std::endl << "  T = " << real_test_statistic << std::endl << std::endl;
      
      // Open the output files
      std::string ti_output_filename_1;
      std::string ti_output_filename_2;
     
      if (ti_output_fn_1) {
        ti_output_filename_1 = args::get(ti_output_fn_1);
      } else {
        ti_output_filename_1 = "Tis.dataset1_" + std::to_string(tis_1.size()) + ".txt";
      }
      if (ti_output_fn_2) {
        ti_output_filename_2 = args::get(ti_output_fn_2);
      } else {
        ti_output_filename_2 = "Tis.dataset2_" + std::to_string(tis_2.size()) + ".txt";
      }
      
      std::ofstream ti_output_file_1;
      ti_output_file_1.open(ti_output_filename_1);
      std::ofstream ti_output_file_2;
      ti_output_file_2.open(ti_output_filename_2);
      
      // Write to output files (both data and Ti values)
      std::cout << "Writing Ti values for dataset 1 to " << ti_output_filename_1 << std::endl;
      for (size_t i=0; i < tis_1.size(); ++i) {
        ti_output_file_1 << std::setprecision(9) << dataset_1[i].s12 << " " << dataset_1[i].s13 << " " << dataset_1[i].s24 << " " << dataset_1[i].s34 << " " << dataset_1[i].s134
                            << " " << tis_1[i] << std::endl;
      }
      std::cout << "Writing Ti values for dataset 2 to " << ti_output_filename_2 << std::endl;
      for (size_t i=0; i < tis_2.size(); ++i) {
        ti_output_file_2 << std::setprecision(9) << dataset_2[i].s12 << " " << dataset_2[i].s13 << " " << dataset_2[i].s24 << " " << dataset_2[i].s34 << " " << dataset_2[i].s134
                            << " " << tis_2[i] << std::endl;
      }
      
      ti_output_file_1.close();
      ti_output_file_2.close();
    }
    // Compute only the test statistic for the current dataset
    else {
      std::cout << "Calculating test statistic for nominal dataset:" << std::endl;
      real_test_statistic = compute_statistic(dataset_1, dataset_2, true);
      std::cout << "  T = " << real_test_statistic << std::endl;
    }
  }
  std::cout << std::endl;
  if (n_permutations) {
    
    // Merge the vectors of events so we can shuffle them
    std::vector<Event> all_events;
    all_events.insert(all_events.end(), dataset_1.begin(), dataset_1.end());
    all_events.insert(all_events.end(), dataset_2.begin(), dataset_2.end());

    size_t N = args::get(n_permutations);

    // Scale the number of events used in permutations
    size_t n_events_1 = dataset_1.size();
    size_t n_events_2 = dataset_2.size();
    if (max_permutation_events_1) {
      n_events_1 = std::min(args::get(max_permutation_events_1), n_events_1);
      n_events_2 = std::round(n_events_1 * ((double) dataset_2.size()/ (double) dataset_1.size()));
    }

    double factor =  1.0*(n_events_1+n_events_2)/(1.0*(dataset_1.size()+dataset_2.size()));
    // Set up the random number generator
    int random_seed = std::mt19937::default_seed;
    if (seed)
      random_seed = args::get(seed);
    std::mt19937 random_generator(random_seed);

    // Open the output files
    
    // File for permuted T values
    std::string permutation_output_filename;
    if (permutation_output_fn) {
      permutation_output_filename = args::get(permutation_output_fn);
    } else {
      permutation_output_filename = "Ts." + std::to_string(dataset_1.size()) + "_";
      permutation_output_filename += std::to_string(dataset_2.size()) + "_";
      permutation_output_filename += std::to_string(N) + "_";
      permutation_output_filename += std::to_string(random_seed) + ".txt";
    }
    std::ofstream output_file;
    if(output_write || permutation_output_fn){
      output_file.open(permutation_output_filename);
    }
    
    // File for min/max permuted Ti values
    std::string permutation_ti_minmax_output_filename;
    if (permutation_ti_minmax_output_fn) {
      permutation_ti_minmax_output_filename = args::get(permutation_ti_minmax_output_fn);
    } else {
      permutation_ti_minmax_output_filename = "Tis_min_max." + std::to_string(dataset_1.size()) + "_";
      permutation_ti_minmax_output_filename += std::to_string(dataset_2.size()) + "_";
      permutation_ti_minmax_output_filename += std::to_string(N) + "_";
      permutation_ti_minmax_output_filename += std::to_string(random_seed) + ".txt";
    }
    std::ofstream ti_minmax_output_file;
    if(output_write || permutation_ti_minmax_output_fn){
      ti_minmax_output_file.open(permutation_ti_minmax_output_filename);
    }
    
    // File for p values
    std::ofstream p_output_file;
    if(!permutations_only){
      p_output_file.open("pvalues.txt", std::iostream::out | std::iostream::app );
    }
    
    int nsig = 0;

    std::cout << "Running " << N << " permutations of " << n_events_1 << " and "
              << n_events_2 << " events using seed " << random_seed << std::endl;
    if(output_write){
      std::cout << "Output filename for permuted T values is " << permutation_output_filename << std::endl;
      if (calculate_ti) std::cout << "Output filename for min/max Ti values is " << permutation_ti_minmax_output_filename << std::endl;
    }

    // Counter to avoid shuffling every time, start large to do a first shuffle
    size_t events_to_skip = all_events.size()+1;

    for (size_t i = 0; i < N; ++i) {
      if ((i+1) % std::max((size_t) 100, N/100) == 0)
        std::cout << "Calculating permutation " << i+1 << " of " << N << std::endl;

      // Reshuffle if we've ran out of events
      if (events_to_skip + n_events_1 + n_events_2 > all_events.size()) {
        std::shuffle(all_events.begin(), all_events.end(), random_generator);
        events_to_skip = 0;
      }

      const std::vector<Event> data_1(all_events.begin()+events_to_skip, all_events.begin()+events_to_skip+n_events_1);
      events_to_skip += n_events_1;
      const std::vector<Event> data_2(all_events.begin()+events_to_skip, all_events.begin()+events_to_skip+n_events_2);
      events_to_skip += n_events_2;
      
      // Calculate permuted T value (and Ti values if specified)
      double test_statistic = -999;
      // Via Ti values
      if (calculate_ti) {
        // Calculate Ti values
        const std::vector<double> perm_tis_1 = compute_statistic_contributions(data_1, data_2);
        const std::vector<double> perm_tis_2 = compute_statistic_contributions(data_2, data_1);
        
        // Sum Ti to get permuted test statistic for this sub-sample
        double total = 0;
        for (size_t i=0; i < perm_tis_1.size(); ++i) total += perm_tis_1[i];
        for (size_t i=0; i < perm_tis_2.size(); ++i) total += perm_tis_2[i];
        test_statistic = total;
        
        // Find the min & max Ti values
        double ti_max = -999;
        double ti_min = 999;
        for (size_t i=0; i < perm_tis_1.size(); ++i) {
          if (perm_tis_1[i] > ti_max) ti_max = perm_tis_1[i];
          if (perm_tis_1[i] < ti_min) ti_min = perm_tis_1[i];
        }
        for (size_t i=0; i < perm_tis_2.size(); ++i) {
          if (perm_tis_2[i] > ti_max) ti_max = perm_tis_2[i];
          if (perm_tis_2[i] < ti_min) ti_min = perm_tis_2[i];
        }
        
        // Write max/min Ti values to file
        if (output_write || permutation_ti_minmax_output_fn) { 
          ti_minmax_output_file << ti_min << " " << ti_max << std::endl;
        }
      }
      // Directly (without Ti values)
      else {
        test_statistic = compute_statistic(data_1, data_2);
      }
      
      // Write permuted T values to file
      if(output_write || permutation_output_fn){output_file << test_statistic << std::endl;}
      
      // Increment counter used for calculating p-value
      if(!permutations_only){
        if (factor * test_statistic > real_test_statistic){nsig++;} 
      }
    }
    if(!permutations_only){
      p_output_file << delta << " "<< (1.0*nsig)/(1.0*N)  << std::endl;
      std::cout<<"p value is: "<<(1.0*nsig)/(1.0*N)  << std::endl;
    }
    p_output_file.close();
    output_file.close();
    ti_minmax_output_file.close();
  }
  
  return 0;
}

int main(int argc, char *argv[]) {
  try {
    return run_energy_test(argc, argv);
  } catch (std::runtime_error& e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }
}
