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

#ifndef LINEAR_FUNCTION_H
#define LINEAR_FUNCTION_H

#include "common.h"
#include "function.h"

namespace apn {
	template <size_t BITS_IN, size_t BITS_OUT>
	class linear_function : public virtual function_in_out_t {
		protected:
			using function_in_out_t::domain;
			using function_in_out_t::codomain;
			size_t domain_dimension = 0;
			size_t codomain_dimension = 0;

		public:
			// using function_in_out_t::function;
			linear_function() {set(0,0);} // linearity constraint
			linear_function(word_out_t y) = delete; // makes no sense unless x == 0

			using function_in_out_t::set;
			using function_in_out_t::test;
			using function_in_out_t::test_codomain;
			using function_in_out_t::get_domain;
			using function_in_out_t::get_codomain;
			size_t get_domain_dimension() const;
			size_t get_codomain_dimension() const;

			void reset();

			// should overload set to add the span
			bitset_in_t add_linear_span(word_in_t x, word_out_t y);
			bitset_in_t add_refining_linear_span(word_in_t x, word_out_t y);
	};
}

using namespace std;
using namespace apn;

template <size_t BITS_IN, size_t BITS_OUT>
void linear_function_in_out_t::reset() {
	function_in_out_t::reset();
	domain_dimension = 0;
	codomain_dimension = 0;
}

template <size_t BITS_IN, size_t BITS_OUT>
size_t linear_function_in_out_t::get_domain_dimension() const {
	return domain_dimension;
}

template <size_t BITS_IN, size_t BITS_OUT>
size_t linear_function_in_out_t::get_codomain_dimension() const {
	return codomain_dimension;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_in_t linear_function_in_out_t::add_linear_span(word_in_t x, word_out_t y) {
	const linear_function_in_out_t& A = *this;
	bitset_in_t span;
	const bitset_in_t domain = A.get_domain();

	domain_dimension++;
	if (!test_codomain(y))
		codomain_dimension++;

	{
		set(x, y);
		span.set(x);
	}
	for (word_in_t pos = 1; ; ++pos) {
		if (domain.test(pos)) {
			word_out_t z = y^A[pos];
			set(x^pos, z);
			span.set(x^pos);
		}

		if (pos == constants_in_t::WORD_MAX)
			break;
	}

	return span;
}

/* Asserts the domain is exactly the set of elements strictly smaller than x (and so x is a power of two) */
template <size_t BITS_IN, size_t BITS_OUT>
bitset_in_t linear_function_in_out_t::add_refining_linear_span(word_in_t x, word_out_t y) {
	const linear_function_in_out_t& A = *this;
	bitset_in_t span = domain << x;

	domain_dimension++;
	if (!test_codomain(y))
		codomain_dimension++;

	{
		set(x, y);
	}
	for (word_in_t pos = 1; pos < x; ++pos) {
		word_out_t z = y^A[pos];
		set(x^pos, z);
	}

	return span;
}

#endif

