// Decider benchmark for the C++ implementation of Bringmann, Künnemann and Nusser,
// "Walking the Dog Fast in Practice" (SoCG 2019),
// https://gitlab.mpi-klsb.mpg.de/anusser/frechet_distance
//
// Usage: decider_bench <curve_dir> <answers_out> <runs> <query_file>...
//
// Times FrechetLight::lessThanWithFilters exactly as deciderComparisonExp in their
// src/paper_experiments.cpp does: the curves of a query file are read first, and only
// the decision calls are timed. Each file is run once untimed (warm-up), then timed
// <runs> times; the minimum is reported. Prints one line per query file:
// "<file> <queries> <yes> <min_ms> <load_ms>", and
// writes "<curve1> <curve2> <distance> <0|1>" per query to <answers_out> (1 when the
// distance of the curves is below the query distance).
#include "frechet_light.h"
#include "parser.h"

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct DeciderQuery
{
	std::string name1, name2;
	Curve curve1, curve2;
	distance_t distance;
};

static std::vector<DeciderQuery> load(std::string const& file, std::string const& dir)
{
	std::ifstream f(file);
	if (!f.is_open()) {
		std::cerr << "cannot open " << file << "\n";
		std::exit(1);
	}
	std::vector<DeciderQuery> queries;
	std::string a, b;
	distance_t d;
	while (f >> a >> b >> d) {
		queries.push_back({a, b, parser::readCurve(dir + a), parser::readCurve(dir + b), d});
	}
	return queries;
}

int main(int argc, char** argv)
{
	if (argc < 5) {
		std::cerr << "usage: decider_bench <curve_dir> <answers_out> <runs> <query_file>...\n";
		return 1;
	}
	std::string dir = argv[1];
	std::ofstream answers(argv[2]);
	int runs = std::atoi(argv[3]);
	answers << std::setprecision(std::numeric_limits<distance_t>::max_digits10);
	std::cout << std::fixed << std::setprecision(3);

	FrechetLight frechet;
	for (int i = 4; i < argc; ++i) {
		auto t0 = std::chrono::high_resolution_clock::now();
		auto queries = load(argv[i], dir);
		double load_ms = std::chrono::duration<double, std::milli>(
			std::chrono::high_resolution_clock::now() - t0).count();
		// Warm-up: one untimed pass over the queries.
		for (auto const& q: queries) { frechet.lessThanWithFilters(q.distance, q.curve1, q.curve2); }
		std::vector<char> result(queries.size());
		double best = std::numeric_limits<double>::infinity();
		for (int run = 0; run < runs; ++run) {
			auto start = std::chrono::high_resolution_clock::now();
			for (std::size_t q = 0; q < queries.size(); ++q) {
				result[q] = frechet.lessThanWithFilters(
					queries[q].distance, queries[q].curve1, queries[q].curve2);
			}
			auto stop = std::chrono::high_resolution_clock::now();
			double ms = std::chrono::duration<double, std::milli>(stop - start).count();
			if (ms < best) { best = ms; }
		}
		std::size_t yes = 0;
		for (std::size_t q = 0; q < queries.size(); ++q) {
			yes += result[q];
			answers << queries[q].name1 << " " << queries[q].name2 << " "
				<< queries[q].distance << " " << int(result[q]) << "\n";
		}
		std::cout << argv[i] << " " << queries.size() << " " << yes << " " << best << " " << load_ms << "\n";
	}
	return 0;
}
