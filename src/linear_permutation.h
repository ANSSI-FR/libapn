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

#ifndef LINEAR_PERMUTATION_H
#define LINEAR_PERMUTATION_H

#include "common.h"
#include "linear_function.h"
#include "permutation.h"

namespace apn {

	template <size_t BITS>
	class linear_permutation : public permutation_t, public linear_function_t {
		protected:
			using linear_function_t::domain;
			using linear_function_t::codomain;
			using linear_function_t::domain_dimension;
			using linear_function_t::codomain_dimension;

			using permutation_t::inverse;

		private:
			static xor_t xor_values[1<<BITS];
			static size_t xor_loops[BITS];
			static size_t xor_words[BITS];

		public:
			static void init_xor_tables();
			using array_t::data;

			linear_permutation() {set(0,0);} // linearity constraint, uses virtual function

			using permutation_t::set;
			using linear_function_t::reset;

			using permutation_t::test;
			using permutation_t::test_codomain;

			using linear_function_t::get_domain;
			using linear_function_t::get_codomain;

			bitset_t add_refining_linear_span(word_t x, word_t y);
			bitset_t quickly_add_refining_linear_span(word_t x, word_t y, bool compute_inverse = true);
			linear_permutation_t compute_inverse() const;
	};
}

using namespace std;
using namespace apn;

template <size_t BITS>
void linear_permutation_t::init_xor_tables() {
	size_t i = 0;
	for (size_t j = 1; j < constants_t::XOR_BITS/constants_t::WORD_BITS; j <<= 1) {
		xor_words[i] = j;
		i++;
	}
	for (; i < BITS; i++) {
		xor_words[i] = constants_t::XOR_BITS/constants_t::WORD_BITS;
	}
	for (i = 0; i < BITS; i++) {
		xor_loops[i] = ((1<<i) / xor_words[i]);
	}	
	for (word_t x = 0; ; x++) {
		xor_t y = 0;
		for (i = 0; i < xor_words[BITS-1]; i++) {
			y <<= constants_t::WORD_BITS;
			y ^= x;
		}
		xor_values[x] = y;
		if (x == constants_t::WORD_MAX)
			break;
	}
}

template <size_t BITS>
linear_permutation_t linear_permutation_t::compute_inverse() const {
	linear_permutation_t inv;
	for (word_t x = 0; x < (1 << domain_dimension); x++)
		inv[(*this)[x]] = x;
	inv.set_domain(get_codomain());
	inv.set_codomain(get_domain());
	return inv;
}

template <size_t BITS>
bitset_t linear_permutation_t::add_refining_linear_span(word_t x, word_t y) {
	const linear_permutation_t& A = *this;
	bitset_t span = domain << x;

	domain_dimension++;
	codomain_dimension++;

	{
		set(x, y);
	}
	for (word_t pos = 1; pos < x; ++pos) {
		word_t z = y^A[pos];
		set(x^pos, z);
	}

	return span;
}

template <size_t BITS>
bitset_t linear_permutation_t::quickly_add_refining_linear_span(word_t x, word_t y, bool compute_inverse) {
	bitset_t span = domain << x;

	xor_t* zptr = (xor_t*) data();
	xor_t* xptr = (xor_t*) (data() + x);
	xor_t yxor = xor_values[y];

	for (size_t i = 0; i < xor_loops[domain_dimension]; i++) {
		//*xptr == _mm_xor_si128(*zptr, yxor);
		*xptr = (*zptr) ^ yxor;
		// Better to do?
		for (size_t j = 0; j < xor_words[domain_dimension]; j++) {
			word_t y = *(((word_t*) xptr) + j);
			codomain.set(y);
			/* Delay computation of inverse? */ 
			if (compute_inverse)
				inverse[y] = x+i*xor_words[domain_dimension]+j;
		}
		xptr++;
		zptr++;
	}

	domain |= span;
	domain_dimension++;
	codomain_dimension++;

	return span;
}

#endif

