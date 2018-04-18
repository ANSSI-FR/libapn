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

#ifndef _INVOLUTION_H_
#define _INVOLUTION_H_

#include <initializer_list>
#include <vector>

#include "apn.h"
#include "common.h"

namespace apn {

	template <size_t BITS>
	class linear_permutation;

	template <size_t BITS>
	class involution : public permutation_t {
		protected:
			using function_t::domain;
			using function_t::codomain;

		public:
			involution() = default; // default so not inherited
			involution(const array_t& a); // calls virtual function set
			involution(std::initializer_list<word_t> list); // calls virtual function set

			void set(const word_t x, const word_t y);
			void reset(word_t x);

			using function_t::get_domain;
			using function_t::get_codomain;

			using function_t::test;
			using function_t::test_codomain;

			bool is_weakly_linearly_minimal() const;
			bool is_weakly_linearly_smaller(const involution_t& minS) const;

			bool minimally_extend_representative(
				const involution_t& minS,
				linear_permutation_t& A,
				bitset_t& nA,
				involution_t& R,
				linear_permutation_t* minA = NULL,
				involution_t* minR = NULL,
				bool early_abort = true
			) const;

			bool propagate_na_constraints(
				const involution_t& minS,
				linear_permutation_t& A,
				bitset_t& nA,
				involution_t& R,
				linear_permutation_t* minA = NULL,
				involution_t* minR = NULL,
				bool early_abort = true
			) const;


			bool make_guess(
				const involution_t& minS,
				linear_permutation_t& A,
				bitset_t& nA,
				involution_t& R,
				linear_permutation_t* minA = NULL,
				involution_t* minR = NULL,
				bool early_abort = true
			) const;

			bool complete_representative(
				const involution_t& minS,
				linear_permutation_t& A,
				involution_t& R,
				linear_permutation_t* minA = NULL,
				involution_t* minR = NULL,
				bool early_abort = true
			) const;


	};
}

using namespace std;
using namespace apn;

template <size_t BITS>
involution_t::involution(initializer_list<word_t> list) {
	const word_t* it = begin(list);
	for (word_t x = 0; x < list.size(); ++x) {
		set(x, *it);
		++it;
		if (x == constants_t::WORD_MAX)
			break;
	}
}

template <size_t BITS>
involution_t::involution(const array_t& a) {
	for (word_t x = 0; ; ++x) {
		set(x, a[x]);
		if (x == constants_t::WORD_MAX)
			break;
	}
}

template <size_t BITS>
void involution_t::set(const word_t x, const word_t y) {
	permutation_t::set(x, y);
	permutation_t::set(y, x);
}

template <size_t BITS>
void involution_t::reset(word_t x) {
	if (test(x)) {
		word_t y = (*this)[x];
		domain.reset(x);
		domain.reset(y);
		codomain.reset(x);
		codomain.reset(y);
	}
}

/* Return value is for early abortion */
template <size_t BITS>
bool involution_t::make_guess(
	const involution_t& minS,
	linear_permutation_t& A,
	bitset_t& nA,
	involution_t& R,
	linear_permutation_t* minA,
	involution_t* minR,
	bool early_abort
) const {

	const involution_t& S = *this;

	linear_permutation_t bak_A = A;
	involution_t bak_R = R;

	word_t dim;
	word_t x;

	/* Test all possible values, zero is not a valid value and will break stuff */
	for(word_t guess=1; ; ++guess) {

		/* Nothing to say if S is not defined there */
		if (!S.test(guess)) {
			if (guess == constants_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* A is not a permutation */
		if (A.test_codomain(guess)) {
			if (guess == constants_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* Pick the smallest unassigned element for A */
		dim = A.get_domain_dimension();
		x = 1 << dim;

		/* Make the guess */
		nA = A.quickly_add_refining_linear_span(x, guess, true);

#if DEBUG
		cout << "--------------------------------------\n";
		cout << "---- extending A ---------------------\n";
		cout << "--------------------------------------\n";
		cout << "S = " << S << endl;
		cout << "minS = " << minS << endl;
		cout << "A = " << A << endl;
		cout << "R = " << R << endl;
		cout << "x = " << (long) x << endl;
		cout << "guess = " << (long) guess << endl;
#endif
		/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
		if (minimally_extend_representative(minS, A, nA, R, minA, minR, early_abort))
			return true;

		if (guess == constants_t::WORD_MAX)
			break;

		/* Then try next guess */
		A = bak_A;
		R = bak_R;

	} /* Loop on guesses */

	return false;
}

template <size_t BITS>
bool involution_t::minimally_extend_representative(
	const involution_t& minS,
	linear_permutation_t& A,
	bitset_t& nA,
	involution_t& R,
	linear_permutation_t* minA,
	involution_t* minR,
	bool early_abort
) const {

	const involution_t& S = *this;

	/* The caller can have set something in nA, treat it */
	if (nA.any()) {
		return propagate_na_constraints(minS, A, nA, R, minA, minR, early_abort);
	}

	/* Are finished? */
	if (((~A.get_domain()) & S.get_domain()).none()) {
		return complete_representative(minS, A, R, minA, minR, early_abort);
	}

	return make_guess(minS, A, nA, R, minA, minR, early_abort);
}

template <size_t BITS>
bool involution_t::propagate_na_constraints(
	const involution_t& minS,
	linear_permutation_t& A,
	bitset_t& nA,
	involution_t& R,
	linear_permutation_t* minA,
	involution_t* minR,
	bool early_abort
) const {
	const involution_t& S = *this;

	word_t x;
	word_t y;
	word_t z;
	size_t dim;

	/* A part */
	while (nA.any()) {
		x = ffs(nA);
		y = A[x];

		/* If S is not defined at A[x], nothing to say */
		if (!S.test(y)) {
			nA.reset(x);
			continue;
		}
	
		z = S[y];

		if (A.test_codomain(z)) {

			y = A.invert(z);
			nA.reset(x);

			if (R.test_codomain(y))
				continue;
			
			R.set(x, y);
			if (!R.is_involution())
				continue;

		} else {

			/* Next undefined value for A */
			dim = A.get_domain_dimension();
			y = (1 << dim);

			if (R.test(y))
				continue;

			/* Update nA */
			nA.reset(x);

			linear_permutation_t tmp = A;
			bitset_t span = tmp.quickly_add_refining_linear_span(y, z, true);

			if (!tmp.is_involution())
				continue;

			/* R is an involution */
			R.set(x, y);
			// R.set(y, x);

			A = tmp;
			nA |= span;

		}

		/* Early abort? */
		if (early_abort && R < minS) {
			return true;
		}

		/* No need to recurse if the result is larger than the current minimum */
		if (R > minS || (minR != NULL && R > *minR)) {
			return false;
		}
	}

	/* Now nA is empty and nB surely is not. */
	return minimally_extend_representative(minS, A, nA, R, minA, minR, early_abort);
}

template <size_t BITS>
bool involution_t::complete_representative(
	const involution_t& minS,
	linear_permutation_t& A,
	involution_t& R,
	linear_permutation_t* minA,
	involution_t* minR,
	bool early_abort
) const {
	const involution_t& S = *this;

	/* Complete R */
	word_t dim = A.get_domain_dimension();
	for (word_t x = 0; x<(1<<dim); ++x) {
		word_t y = A[x];
		if (!R.test(x) && S.test(y)) {
			word_t z = S[y];
			if (A.test_codomain(z)) {
				R.set(x, A.invert(z));
			}
		}
	}/*x*/

	if (early_abort && R < minS)
		return true;

	/* Did we find a smaller representative? */
	if (minR != NULL && !(*minR < R)) {
		*minR = R;
		*minA = A;
	}

	return false;
}

#endif


template <size_t BITS>
bool involution_t::is_weakly_linearly_minimal() const {
	const involution_t& S = *this;

	if (test(0)) {
		if (S[0] != 0) {
			if (S[0] != 1)
				return false;

			/* Here, S[0]=1 */
			if (S[1] != 0)
				return false; /* since S is an involution */

		} else {
			if (test(1) && S[1] > 1)
				return false;
		}
	}

	return !is_weakly_linearly_smaller(S);
}


template <size_t BITS>
bool involution_t::is_weakly_linearly_smaller(const involution_t& minS) const {

	involution_t R;

    /* A is a linear involution */
    linear_permutation_t A;

    bitset_t nA = 0;
    nA.set(0);

    // if (S[0]==0) {
    //     R.set(0, 0);

    // } else {
    //     // R.set(0, 1);
    //     // R.set(1, 0);

    //     // A.set(1, S[0]);
    //     // A.set(S[0], 1);

    //     // nA.set(0);
    //     // nA.set(1);
    // }

	if (minimally_extend_representative(minS, A, nA, R, NULL, NULL, true))
		return true;

	return false;
}


