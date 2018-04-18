/*
Copyright 2017 Jean-Pierre Flori, Jérémy Jean

This file is part of libapn.

libapn is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

libapn is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with libapn; if not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <chrono>
#include <ctime>
#include <unistd.h>

#include <queue>

#include <pthread.h>

#include "config.h"
#include "apn.h"

#define BITS_IN 6
#define BITS_OUT 4
#define DIFF_UNIFORMITY 6
#define DEGREE 3

#define NUM_THREADS 24
#define SHARE_THRESHOLD 1

using namespace apn;
using namespace std;

static chrono::system_clock::time_point start;

static queue<function_in_out_t> sboxes;
static pthread_mutex_t mutex_sboxes = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond_sboxes = PTHREAD_COND_INITIALIZER;
static int wanted_sboxes = 0;
//static pthread_mutex_t mutex_wanted = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond_wanted = PTHREAD_COND_INITIALIZER;

static pthread_mutex_t mutex_print = PTHREAD_MUTEX_INITIALIZER;

template<>
xor_in_t apn::linear_permutation_in_t::xor_values[1<<BITS_IN] = {};
template<>
size_t apn::linear_permutation_in_t::xor_loops[BITS_IN] = {};
template<>
size_t apn::linear_permutation_in_t::xor_words[BITS_IN] = {};

template<>
xor_out_t apn::linear_permutation_out_t::xor_values[1<<BITS_OUT] = {};
template<>
size_t apn::linear_permutation_out_t::xor_loops[BITS_OUT] = {};
template<>
size_t apn::linear_permutation_out_t::xor_words[BITS_OUT] = {};

/*
** Recursively add guesses to the partially defined Sbox 
** and early-abort Sbox prefixes depending on:
**		- bijectivity
**		- APN-ness
**		- Affine equivalence 
*/
void traverse_functions_rec(function_in_out_t& S) {

	word_in_t position = S.get_position();

	/* Check whether that position has already been determined */
	if (S.test(position)) {

		/* Display progress of the computation */
		if (position <= 8) { // || (position == 24) || (position == 30)) {
			//pthread_mutex_lock(&mutex_print);
			cerr << "[-] Exploring " << S << " at " << &S << " by " << pthread_self() << " pid " << getpid() << endl;
			//cerr.flush();
			//pthread_mutex_unlock(&mutex_print);
		}

		/* (APN)
		 * Check that the newly added value does not contract the APN property
		 * Warning: ddt is modified by checkAPN
		 */
#if (DIFF_UNIFORMITY==2)
	    if (!S.extend_bbt(true)) {
			return ;
		}

#else /*(DIFF_UNIFORMITY>2)*/
		if (!S.extend_ddt(true)) {
			//cerr << "wrong ddt" << endl;
			return ;
		}
#endif

		/* (Affine equivalence)
         * Only continue with the partial permutation S if there is one completion
         * of it that is the representant of the affine equivalence class
         */
		if ((position < (1<<BITS_IN)/2 || position == constants_in_t::WORD_MAX) && !S.is_affinely_minimal()) {
		// if (!S.is_affinely_minimal()) {
			//cerr << "not minimal" << endl;
			return;
		}

		//cerr << "Passed tests" << endl;
		/* If finished */
		if (position == constants_in_t::WORD_MAX) {
			/* Here, we have a leaf in the tree that passes all the criteria */
			pthread_mutex_lock(&mutex_print);
			chrono::system_clock::time_point end = chrono::system_clock::now();
			time_t end_time = chrono::system_clock::to_time_t(end);
			chrono::duration<double> elapsed = end-start;
			cout << "[+] Found function: [max=" << (int) S.compute_ddt() << ", deg=" << (int) S.compute_degree() << ", inv=" << (int) S.is_involution() << "] ";
			cout << S.to_array_string() << endl;
			S.print_ddt();
			cout << "[#] At " << ctime(&end_time);
			cout << "[#] Time elapsed: " << elapsed.count() << " seconds" << endl;
			cout.flush();
			pthread_mutex_unlock(&mutex_print);
			return;
		}

		S.inc_position();
		traverse_functions_rec(S);
		return;
	}


	/*
	 * If we have a filter on the degree of the permutation,
	 * check whether a constraint of the degree can yield S[position]
	 */
#if DEGREE != 0
	/* Give the last set position to check whether we can directly determine more positions */
	if (S.propagate_degree_constraints(position)) {
		traverse_functions_rec(S);
		return;
	}
#endif

	if (__builtin_popcount(position) == 1) {

		/* Enforced by EA equivalence */
		S.set(position, 0);
		traverse_functions_rec(S);

	} else {

		word_out_t guess = 0;

		if (wanted_sboxes) {
			
			/* wanted_sboxes could have changes but that is ok */
			pthread_mutex_lock(&mutex_sboxes);
			for (; guess < constants_out_t::WORD_MAX-SHARE_THRESHOLD; guess++) {
				function_in_out_t recS = S;

				recS.set(position, guess);
				// cerr << pthread_self() << " giving " << recS << endl;
				sboxes.push(recS);
			}
			function_in_out_t recS = S;
			recS.set(position, guess);
			// cerr << pthread_self() << " keeping " << recS << endl;
			pthread_cond_broadcast(&cond_sboxes);
			pthread_mutex_unlock(&mutex_sboxes);
		}

		/* At least one left for us */
		for (; ; ++guess) {
			function_in_out_t recS = S;

			/* Perform the guess */
			recS.set(position, guess);
			traverse_functions_rec(recS);

			if (guess == constants_out_t::WORD_MAX)
				break;
		}/*guess*/
	}
}/*traverse_functions_rec*/


static void* traverse_functions(void*) {
	function_in_out_t S;
	bool waited = false;
	while (true) {
		pthread_mutex_lock(&mutex_sboxes);
		if (sboxes.empty()) {
			//pthread_mutex_lock(&mutex_wanted);
			//cerr << pthread_self() << " waiting" << endl;
			wanted_sboxes++;
			waited = true;
			pthread_cond_signal(&cond_wanted);
			//pthread_mutex_unlock(&mutex_wanted);
			pthread_cond_wait(&cond_sboxes, &mutex_sboxes);
		}
		while (sboxes.empty()) {
			//cerr << pthread_self() << " waiting" << endl;
			pthread_cond_wait(&cond_sboxes, &mutex_sboxes);
		}
		S = sboxes.front();
		//cerr << pthread_self() << " taking " << S << endl;
		sboxes.pop();
		//S = sboxes.back();
		//sboxes.pop_back();
		//pthread_mutex_lock(&mutex_wanted);
		if (waited) {
			wanted_sboxes--;
			waited = false;
		}
		//pthread_mutex_unlock(&mutex_wanted);
		pthread_mutex_unlock(&mutex_sboxes);
		traverse_functions_rec(S);
	}
}

/*
 * Perform a depth-first tree-based search of all BITS-bit APN functions,
 * while pruning branches (with affine equivalence and APN checks).
 * Note: the algorithm does not consider general functions, but only bijective ones.
 */
int main(void) {

	/* Log parameters of the search */
	cout << "[-] Parameters:" << endl;
	cout << "         NUM_THREADS: " << NUM_THREADS      << endl;
	cout << "             BITS_IN: " << BITS_IN          << endl;
	cout << "            BITS_OUT: " << BITS_OUT         << endl;
	cout << "              DEGREE: " << DEGREE           << endl;
	cout << "     DIFF_UNIFORMITY: " << DIFF_UNIFORMITY  << endl;
	cout << endl;
	
	/* Start the timer */
	start = chrono::system_clock::now();
	time_t start_time = chrono::system_clock::to_time_t(start);
	cout << "[#] Starting computation at " << ctime(&start_time);

	linear_permutation_in_t::init_xor_tables();
	linear_permutation_out_t::init_xor_tables();

	function_in_out_t S;
	S.set(0,0);
	S.inc_position();

	/* Search for DIFF_UNIFORMITY-uniform permutation */
	S.set_max_uniformity(DIFF_UNIFORMITY);

	/* Search for functions having at most that degree */
	S.set_max_degree(DEGREE);

	/* Add the Sbox to the container */
	sboxes.push(S);
	//cerr << pthread_self() << " pushing " << S << endl;

	pthread_attr_t attr;
	pthread_t workers[NUM_THREADS];

	pthread_attr_init(&attr);
	size_t stacksize;
	stacksize = 128*1024*1024;
	pthread_attr_setstacksize(&attr, stacksize);
	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_create(&workers[i], &attr, traverse_functions, NULL);
	}
	pthread_attr_destroy(&attr);

	pthread_mutex_lock(&mutex_sboxes);
	while (wanted_sboxes != NUM_THREADS) {
		//pthread_cond_wait(&cond_wanted, &mutex_wanted);
		pthread_cond_wait(&cond_wanted, &mutex_sboxes);
	}
	pthread_mutex_unlock(&mutex_sboxes);

	//if (!sboxes.empty()) abort();
	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_cancel(workers[i]);
	}

	/* Stop the timer */
	chrono::system_clock::time_point end = chrono::system_clock::now();
	time_t end_time = chrono::system_clock::to_time_t(end);
	chrono::duration<double> elapsed = end-start;
	cout << "[#] Finished computation at " << ctime(&end_time);
	cout << "[#] Time elapsed: " << elapsed.count() << " seconds" << endl;

	return 0;
}
