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

#include <queue>

#include <pthread.h>

#include "config.h"
#include "apn.h"

#define BITS 5
#define DIFF_UNIFORMITY 2
#define DEGREE 0

#define NUM_THREADS 24
#define SHARE_THRESHOLD 1

using namespace apn;
using namespace std;

static chrono::system_clock::time_point start;

static queue<function_t> sboxes;
static pthread_mutex_t mutex_sboxes = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond_sboxes = PTHREAD_COND_INITIALIZER;
static int wanted_sboxes = 0;
//static pthread_mutex_t mutex_wanted = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cond_wanted = PTHREAD_COND_INITIALIZER;

static pthread_mutex_t mutex_print = PTHREAD_MUTEX_INITIALIZER;

template<>
xor_t apn::linear_permutation_t::xor_values[1<<BITS] = {};
template<>
size_t apn::linear_permutation_t::xor_loops[BITS] = {};
template<>
size_t apn::linear_permutation_t::xor_words[BITS] = {};

/*
** Recursively add guesses to the partially defined Sbox 
** and early-abort Sbox prefixes depending on:
**		- APN-ness
**		- Affine equivalence 
*/
void traverse_functions_rec(function_t& S) {
	size_t position = S.get_position();

	/* Display progress of the computation */
	if (position <= 7) { // || (position == 24) || (position == 30)) {
		//pthread_mutex_lock(&mutex_print);
		cerr << "[-] Exploring " << S << endl;
		//pthread_mutex_unlock(&mutex_print);
	}

	/* Check whether that position has already been determined */
	if (S.test(position)) {

		/* (APN)
		 * Check that the newly added value does not contract the APN property
		 * Warning: ddt is modified by checkAPN
		 */
#if (DIFF_UNIFORMITY==2)
		if (!S.extend_bbt(true)) return ;

#else /*(DIFF_UNIFORMITY>2)*/
		if (!S.extend_ddt(true)) return ;
#endif

		/* (Affine equivalence)
         * Only continue with the partial function S if there is one completion
         * of it that is the representant of the affine equivalence class
         */
		if ((position < (1<<BITS)/2 || position ==constants_t::WORD_MAX -1) && !S.is_affinely_minimal()) {
			return;
		}

		/* If finished */
		if (position == constants_t::WORD_MAX) {
			if (!S.is_extended_affinely_minimal()) return;
			/* Here, we have a leaf in the tree that passes all the criteria */
			pthread_mutex_lock(&mutex_print);
			chrono::system_clock::time_point end = chrono::system_clock::now();
			time_t end_time = chrono::system_clock::to_time_t(end);
			chrono::duration<double> elapsed = end-start;
			cout << "[+] Found APN: [max=" << (int)S.compute_ddt() << ", deg=" << (int)S.compute_degree() << "] ";
			cout << S.to_array_string() << endl;
			cout << "[#] At " << ctime(&end_time);
			cout << "[#] Time elapsed: " << elapsed.count() << " seconds" << endl;
			//cout << "[#] EA minimal: " << S.is_extended_affinely_minimal() << endl;
			cout.flush();
			pthread_mutex_unlock(&mutex_print);
			return ;
		}

		/* Continue the tree search after that position */
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

	/*
	 * For a power-of-2 position
	 */
	if (__builtin_popcount(position) == 1) {

		/* Enforced by EA equivalence */
		S.set(position, 0);
		traverse_functions_rec(S);

	} else {

		word_t guess = 0;

		if (wanted_sboxes) {
			
			/* wanted_sboxes could have changes but that is ok */
			pthread_mutex_lock(&mutex_sboxes);
			for (; guess < constants_t::WORD_MAX-SHARE_THRESHOLD; guess++) {
				function_t recS = S;

				recS.set(position, guess);
				sboxes.push(recS);
			
				if (guess == constants_t::WORD_MAX)
					break;
			}
			pthread_cond_broadcast(&cond_sboxes);
			pthread_mutex_unlock(&mutex_sboxes);
		}

		/* At least one left for us */
		for (; ; ++guess) {
			function_t recS = S;

			/* Perform the guess */
			recS.set(position, guess);
			traverse_functions_rec(recS);

			if (guess == constants_t::WORD_MAX)
				break;
		}/*guess*/
	}
}


static void* traverse_functions(void*) {
	function_t S;
	bool waited = false;
	while (true) {
		pthread_mutex_lock(&mutex_sboxes);
		if (sboxes.empty()) {
			//pthread_mutex_lock(&mutex_wanted);
			wanted_sboxes++;
			waited = true;
			pthread_cond_signal(&cond_wanted);
			//pthread_mutex_unlock(&mutex_wanted);
			pthread_cond_wait(&cond_sboxes, &mutex_sboxes);
		}
		while (sboxes.empty()) {
			pthread_cond_wait(&cond_sboxes, &mutex_sboxes);
		}
		S = sboxes.front();
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
 */
int main(void) {

	setlinebuf(stdout);
	setlinebuf(stderr);

	/* Log parameters of the search */
	cout << "[-] Parameters:" << endl;
	cout << "         NUM_THREADS: " << NUM_THREADS      << endl;
	cout << "                BITS: " << BITS             << endl;
	cout << "              DEGREE: " << DEGREE           << endl;
	cout << "     DIFF_UNIFORMITY: " << DIFF_UNIFORMITY  << endl;
	cout << endl;
	
	/* Start the timer */
	start = chrono::system_clock::now();
	time_t start_time = chrono::system_clock::to_time_t(start);
	cout << "[#] Starting computation at " << ctime(&start_time);

	linear_permutation_t::init_xor_tables();

	function_t S;
	S.set(0,0);
	S.inc_position();

	/* Search for DIFF_UNIFORMITY-uniform sboxes */
	S.set_max_uniformity(DIFF_UNIFORMITY);

	/* Search for permutations having at most that degree */
	S.set_max_degree(DEGREE);

	sboxes.push(S);

	pthread_t workers[NUM_THREADS];

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_create(&workers[i], NULL, traverse_functions, NULL);
	}

	while (wanted_sboxes != NUM_THREADS) {
		//pthread_cond_wait(&cond_wanted, &mutex_wanted);
		pthread_cond_wait(&cond_wanted, &mutex_sboxes);
	}

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
