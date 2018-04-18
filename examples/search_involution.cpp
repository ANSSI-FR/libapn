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

#define BITS 6
#define DIFF_UNIFORMITY 2
#define DEGREE 4

#define NUM_THREADS 24
#define SHARE_THRESHOLD 4

using namespace apn;
using namespace std;

static chrono::system_clock::time_point start;

static queue<involution_t> sboxes;
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
**		- bijectivity
**		- APN-ness
**		- Affine equivalence 
*/
void traverse_involutions_rec(involution_t& S) {
	size_t position = S.get_position();

	/* Display progress of the computation */
	if (position <= 10) { // || (position == 24) || (position == 30)) {
		// pthread_mutex_lock(&mutex_print);
		cerr << "[-] Exploring " << S << endl;
		// pthread_mutex_unlock(&mutex_print);
	}

	/* Check whether that position has already been determined */
	if (S.test(position)) {

		/* (APN)
		 * Check that the newly added value does not contradict the APN property
		 * Warning: ddt is modified by checkAPN
		 */
#if (DIFF_UNIFORMITY==2)
		if (!S.extend_bbt(true)) return ;
#else /*(DIFF_UNIFORMITY>2)*/
		if (!S.extend_ddt(true)) return ;
#endif

		if ((position < (1<<BITS)/2 || position == constants_t::WORD_MAX) && !S.is_weakly_linearly_minimal()) return ;
		// if (!S.is_linearly_minimal()) return ;

		/* If finished */
		if(position == constants_t::WORD_MAX) {
			/* Here, we have a leaf in the tree that passes all the criteria */
			pthread_mutex_lock(&mutex_print);
			chrono::system_clock::time_point end = chrono::system_clock::now();
			time_t end_time = chrono::system_clock::to_time_t(end);
			chrono::duration<double> elapsed = end-start;
			cout << "[+] Found APN: [max=" << (int) S.compute_ddt() << ", deg=" << (int) S.compute_degree() << ", inv=" << (int) S.is_involution() << "] ";
			cout << S.to_array_string() << endl;
			cout << "[#] At " << ctime(&end_time);
			cout << "[#] Time elapsed: " << elapsed.count() << " seconds" << endl;
			cout.flush();
			pthread_mutex_unlock(&mutex_print);
	
			return ;
		}

		S.inc_position();
		traverse_involutions_rec(S);
		return;
	}

// #if 0

	/*
	 * If we have a filter on the degree of the permutation,
	 * check whether a constraint of the degree can yield S[position]
	 */
#if DEGREE != 0
	/* Give the last set position to check whether we can directly determine more positions */
	if (S.propagate_degree_constraints(position)) {
		traverse_involutions_rec(S);
		return;
	}
#endif

	vector<word_t> guesses;
	unsigned int i = 0;

	for (word_t guess = position; ; ++guess) {
		/* (Bijectivity) 
		 * Check that "guess" is not present already in S[2..position-1]
		 */
		if(S.test_codomain(guess)) {
			if (guess == constants_t::WORD_MAX)
				break;
			else
				continue;
		}
		guesses.push_back(guess);
		if (guess == constants_t::WORD_MAX)
			break;
	}

	if (wanted_sboxes && guesses.size() > SHARE_THRESHOLD) {
		// wanted_sboxes could have changes but that is ok
		pthread_mutex_lock(&mutex_sboxes);
		for (; i < guesses.size()-SHARE_THRESHOLD; i++) {
			word_t guess = guesses[i];
			involution_t recS = S;

			recS.set(position, guess);

			sboxes.push(recS);
			//sboxes.push_back(recS);
		}
		pthread_cond_broadcast(&cond_sboxes);
		pthread_mutex_unlock(&mutex_sboxes);
	}

	// At least one left for us
	for (; i < guesses.size(); ++i) {
		word_t guess = guesses[i];
		involution_t recS = S;

		/* Perform the guess */
		recS.set(position, guess);

		traverse_involutions_rec(recS);
	}/*guess*/

}/*traverse_permutations_rec*/


static void* traverse_involutions(void*) {
	involution_t S;
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
		traverse_involutions_rec(S);
	}
}

void test() {

	linear_permutation_t::init_xor_tables();


    // involution_t S = {12,10,13,3,14,11,15,7,8,9,1,5,0,2,4,6};
    // involution_t S = {1, 0, 2, 4, 3, 5, 8, 10, 6, 13, 7, 11, 14, 9, 12, 15};

	involution_t S = {0,19,14,61,47,36,49,26,33,58,11,10,56,39,2,15,29,20,40,1,17,41,22,30,62,31,7,60,28,16,23,25,55,8,37,35,5,34,51,13,18,21,52,48,46,45,44,4,43,6,50,38,42,54,53,32,12,57,9,63,27,3,24,59};
    
    // involution_t S = {0,1,13,3,14,11,15,7,8,9,10,5,12,2,4,6};
	// involution_t S = {0, 1, 2, 3, 4, 5, 8, 10, 6, 9, 7, 14, 13, 12, 11, 15};
    
    cout << S.is_weakly_linearly_minimal() << " <-- S.is_linearly_minimal()" << endl;
    cout << S.is_involution() << " <-- S.is_involution()" << endl;
    cout << S.compute_degree() << " <-- S.compute_degree()" << endl;
	cout << S.compute_ddt() << " <-- S.compute_ddt()" << endl;


}



/*
 * Perform a depth-first tree-based search of all BITS-bit APN permutations,
 * while pruning branches (with affine equivalence and APN checks).
 * Note: the algorithm does not consider general functions, but only bijective ones.
 */
int main(void) {

	// test(); return 0;

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

	involution_t S;
	S.set(0,0);
	S.inc_position();

	/* To converge towards Dillon's involution */
	// S.set(1,19); S.set(19,1); S.inc_position();
	// S.set(2,14); S.set(14,2); S.inc_position();
	// S.set(3,61); S.set(61,3); S.inc_position();
	// S.set(4,47); S.set(47,4); S.inc_position();
	// S.set(5,36); S.set(36,5); S.inc_position();
	// S.set(6,49); S.set(49,6); S.inc_position();
	// S.set(7,26); S.set(26,7); S.inc_position();
	// S.set(8,33); S.set(33,8); S.inc_position();
	// S.set(9,58); S.set(58,9); S.inc_position();
	// S.set(10,11); S.set(11,10); S.inc_position(); S.inc_position();
	// S.set(12,56); S.set(56,12); S.inc_position();
	// S.set(13,39); S.set(39,13); S.inc_position(); S.inc_position();
	// S.set(15,15); S.set(15,15); S.inc_position();
	// S.set(16,29); S.set(29,16); S.inc_position();


	/* Search for DIFF_UNIFORMITY-uniform permutation */
	S.set_max_uniformity(DIFF_UNIFORMITY);

	/* Search for involutions having at most that degree */
	S.set_max_degree(DEGREE);

	/* Add the Sbox to the container */	
	sboxes.push(S);

	pthread_t workers[NUM_THREADS];

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_create(&workers[i], NULL, traverse_involutions, NULL);
	}

	pthread_mutex_lock(&mutex_sboxes);
	while (wanted_sboxes != NUM_THREADS) {
		//pthread_cond_wait(&cond_wanted, &mutex_wanted);
		pthread_cond_wait(&cond_wanted, &mutex_sboxes);
	}
	pthread_mutex_unlock(&mutex_sboxes);

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
