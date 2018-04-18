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

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include <initializer_list>
#include <vector>
// #include <map>

#include "common.h"
#include "function.h"

namespace apn {
	template <size_t BITS>
	class linear_permutation;

	template <size_t BITS>
	class permutation : public virtual function_t {
		protected:
			using function_t::domain;
			using function_t::codomain;

			array_t inverse;

			// std::map<
			// 	word_t,				 // The position that can be determined with previous values at smaller positions
			// 	std::vector<word_t>	 // The list of positions to sum
			// > relations; // Relations (vars,pos) that allows to compute S[pos] from the variables in vars

		public:
			permutation() = default; // default so not inherited, calls function default constructor
			permutation(word_t y) = delete; // makes no sense
			permutation(const array_t& a); // has to be copied: calls virtual function set
			permutation(std::initializer_list<word_t> list); // has to be copied: calls virtual function set

			void set(const word_t x, const word_t y);
			using function_t::reset;
			void reset(word_t x);

			using function_t::get_domain;
			using function_t::get_codomain;

			using function_t::test;
			using function_t::test_domain;
			using function_t::test_codomain;

			bitset_t preimages(word_t y) const;
			bitset_t inverse_image(const bitset_t& set) const;
			word_t invert(word_t y) const;
			permutation_t get_inverse() const;

			using function_t::minimally_extend_representative;
			bool propagate_nb_constraints(const function_t& minS,
				linear_permutation_t& A, linear_permutation_t& B,
				bitset_t& nA, bitset_t& nB,
				function_t& R,
				linear_permutation_t* minA, linear_permutation_t* minB,
				function_t* minR,
				bool early_abort = false) const;

			bool is_affinely_smaller(const function_t& minS) const;

			permutation_t get_linear_representative() const;
			permutation_t get_affine_representative() const;

			// std::pair<bool, word_t> propagate_degree_constraints(const word_t position) const;
	};
}

using namespace std;
using namespace apn;

template <size_t BITS>
permutation_t::permutation(initializer_list<word_t> list) {
	const word_t* it = std::begin(list);
	for (word_t x = 0; x < list.size(); ++x) {
		set(x, *it);
		++it;
		if (x == constants_t::WORD_MAX)
			break;
	}
}

template <size_t BITS>
permutation_t::permutation(const array_t& a) {
	for (word_t x = 0; ; ++x) {
		set(x, a[x]);
		if (x == constants_t::WORD_MAX)
			break;
	}
}

template <size_t BITS>
void permutation_t::set(const word_t x, const word_t y) {
	/*
	if (test_codomain(y))
		return false;
	*/
	function_t::set(x, y);
	inverse[y] = x;
}

template <size_t BITS>
void permutation_t::reset(word_t x) {
	if (test(x)) {
		word_t y = (*this)[x];
		domain.reset(x);
		codomain.reset(y);
	}
}

template <size_t BITS>
word_t permutation_t::invert(word_t y) const {
	return inverse[y];
}

template <size_t BITS>
bitset_t permutation_t::preimages(word_t y) const {
	bitset_t preimages;

	if (test_codomain(y))
		preimages.set(inverse[y]);

	return preimages;
}

template <size_t BITS>
bitset_t permutation_t::inverse_image(const bitset_t& set) const {
	bitset_t range = get_codomain() & set;
	bitset_t image;

	for (word_t y = 0; ; ++y) {
		if (range[y]) {
			image.set(inverse[y]);
		}
		if (y == constants_t::WORD_MAX)
			break;
	}/*x*/

	return image;
}

template <size_t BITS>
permutation_t permutation_t::get_inverse() const {
	permutation_t inv(inverse);
	inv.set_domain(get_codomain());
	inv.set_codomain(get_domain());
	return inv;
}

template <size_t BITS>
bool permutation_t::propagate_nb_constraints(
	const function_t& minS,
	linear_permutation_t& A, linear_permutation_t& B,
	bitset_t& nA, bitset_t& nB,
	function_t& R,
	linear_permutation_t* minA, linear_permutation_t* minB,
	function_t* minR,
	bool early_abort
) const {
	const permutation_t& S = *this;

	word_t x;
	word_t y;
	word_t z;
	size_t dim;

#if SKIP_CLEANUP
#else
	bitset_t common;
#endif
	/* B part */
	while (nB.any()) {
		y = ffs(nB);
		z = B[y];

#if SKIP_CLEANUP
		if (S.test_codomain(z) && !A.test_codomain(S.invert(z)))
#else
		if (S.test_codomain(z))
#endif
		{
			break;
		}
		nB.reset(y);
	}

	if (nB.none())
		return minimally_extend_representative(minS, A, B, nA, nB, R, minA, minB, minR, early_abort);

	/* Next undefined value for A*/
	dim = A.get_domain_dimension();
	x = (1 << dim);
	R.set(x, y);

	/* Early abort ?*/
	if (early_abort && R < minS) {
		return true;
	}
	/* No need to recurse if the result is larger than the current minimum */
	if (R > minS || (minR != NULL && R > *minR)) {
		return false;
	}

#if SKIP_CLEANUP
	nA = A.quickly_add_refining_linear_span(x, S.invert(z), false);
#else
	nA = A.get_domain() << x;
	common = S.image(A.add_refining_linear_span(x, S.invert(z)));
	/* S was maybe somehow linear and we have to reset y anyway */
	nB &= ~(B.inverse_image(common));
	/* this could remove points in nA where S is not defined */
	//nA &= ~(A.inverse_image(S.inverse_image(common)));
#endif

	/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
	if (minimally_extend_representative(minS, A, B, nA, nB, R, minA, minB, minR, early_abort))
	//if (propagate_na_constraints(minS, A, B, nA, nB, R, minA, minB, minR, early_abort))
		return true;

	return false;
}

/*
template <size_t BITS>
bool permutation_t::quick_non_minimality_check() {
	if (test(0) && test(1)) {
		word_t x, swap;
		word_t last = S.get_domain().count() - 1;
		while (x = 1; x <= last; x <<= 1) dim++;
		swap = x >> 1;
		for (x = 0; x <= swap; x++) {
			if (!test(x))
				return false;
		}

		if (!test(last))
			return false;
		if (S[swap] < S[last])
			return true;
	}
	return false;
}
*/

template <size_t BITS>
bool permutation_t::is_affinely_smaller(const function_t& minS) const {
	const permutation_t& S = *this;

	if (minS.test(0) && minS[0] != 0)
		return true;

	/* Iterate through shifts of the original S-box with 0 in zero position */
	for (word_t a = 0; ; ++a) {
		if (!test(a)) {
			if (a == constants_t::WORD_MAX)
				break;
			else
				continue;
		}

		permutation_t shifted_function;
		/* Set shifted[0] = 0 */
		word_t b = S[a];

		for (word_t x = 0; ; ++x) {
			word_t y = x^a;
			if (test(y)) {
				shifted_function.set(x, S[y]^b);
			}
			if (x == constants_t::WORD_MAX)
				break;
		}/*x*/

		if (shifted_function.is_linearly_smaller(minS)) {
			return true;
		}

		if (a == constants_t::WORD_MAX)
			break;
	}/*a*/

	return false;
}


/* Assert S[0] is set */
template <size_t BITS>
permutation_t permutation_t::get_linear_representative() const {

	const permutation_t& S = *this;

	linear_permutation_t A;
	linear_permutation_t B;
	permutation_t R;

	linear_permutation_t minA;
	linear_permutation_t minB;
	permutation_t minR;

	bitset_t nA = 0;
	bitset_t nB = 0;

	if (S[0] != 0) {
		nA.set(0);
		nB.set(0); 
	} else {
		R.set(0, 0);
	}

	minimally_extend_representative(S, A, B, nA, nB, R, &minA, &minB, &minR, false);
	return minR;
}


/* Assert S[0] is set */
template <size_t BITS>
permutation_t permutation_t::get_affine_representative() const {

	const permutation_t& S = *this;
	permutation_t minS = S;

	for (word_t a = 0; a<=constants_t::WORD_MAX; ++a) {
		for(word_t b = 0; b<=constants_t::WORD_MAX; ++b) {

			permutation_t shifted_perm;

			for (word_t x = 0; x<=constants_t::WORD_MAX; ++x) {
				shifted_perm.set(x, S[x^a]^b);
			}/*x*/
			
			permutation_t tmp = shifted_perm.get_linear_representative();
			if(tmp<minS) minS = tmp;

		}/*b*/
	}/*a*/

	return minS;
}

#if 0
/*
 * Check whether the values in S[0..position] can be used to derive S[position+1]
 * when the function S is assumed to be of known degree
 */
template <size_t BITS>
std::pair<bool, word_t> permutation_t::propagate_degree_constraints(const word_t position) const {
	const permutation_t& S = *this;

	std::pair<bool, word_t> ret = function_t::propagate_degree_constraints(position);

	/* If S[position+1] cannot be determined using degree relations */
	if (!ret.first) return make_pair(false, 0);

	/*
	 * If the suggested value is already in the image of S, 
	 * then S cannot be a permutation_t of degree "DEGREE"
	 */
	if (test_codomain(ret.second)) return make_pair(false, 0);

	/* Return the value */
	return make_pair(true, ret.second);

}
#endif

#endif


