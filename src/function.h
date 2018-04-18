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

#ifndef FUNCTION_H
#define FUNCTION_H

#include <cstring>

#include <string>
#include <iostream>

#include <initializer_list>
#include <vector>
#include <algorithm>

#include "common.h"

namespace apn {

	template <size_t BITS_IN, size_t BITS_OUT>
	class function;

	template <size_t BITS_IN, size_t BITS_OUT>
	class linear_function;

	template <size_t BITS>
	class permutation;

	template <size_t BITS>
	class linear_permutation;

	template <size_t BITS>
	class involution;

	template <size_t BITS_IN, size_t BITS_OUT>
	class function : public array_in_out_t {
		protected:
			bitset_in_t domain;
			bitset_out_t codomain;
			bbt_t bbt;
			ddt_t ddt = {{}};
			//vector<int> degrees;
			
			bool apn = true;
			int max_uniformity = 0;
			int uniformity = 0;
			int max_degree = 0;
			int degree = 0;
			word_in_t position = 0; // used to left refine functions

		public:
			function() = default;
			function(word_out_t x) {
				fill(x);
				domain.set();
				codomain.set(x);
			}
			function(const array_in_out_t& a);
			function(std::initializer_list<word_out_t> list);
			
			virtual void set(word_in_t x, word_out_t y);
			virtual void reset();
			virtual void reset(word_in_t x);

			bool test(word_in_t x) const;
			bool test_domain(word_in_t x) const;
			bool fix_codomain(word_out_t y);
			bool test_codomain(word_out_t y) const;

			void set_domain(const bitset_in_t set); 
			void set_codomain(const bitset_out_t set); 

			bitset_in_t get_domain() const;
			bitset_out_t get_codomain() const;
			bitset_out_t image(const bitset_in_t set) const;
			virtual bitset_in_t inverse_image(const bitset_out_t& set) const;
			virtual bitset_in_t preimages(word_out_t y) const;
			vector<word_in_t> preimages_vector(word_out_t y) const;

			bool is_permutation() const;
			bool is_involution() const;

			word_in_t get_position() const;
			void set_position(word_in_t);
			void inc_position();

			int get_max_uniformity() const;
			void set_max_uniformity(int uniformity);
			bool extend_bbt(bool early_abort);
			bool extend_ddt(bool early_abort);
			int compute_ddt();
			int get_uniformity() const;
			bool is_apn() const;

			//static table_degree_constraints;
			virtual bool propagate_degree_constraints(word_in_t position);
			void set_max_degree(int degree);
			int get_max_degree() const;
			int compute_degree();
			//std::vector<int> get_degrees() const;
			int get_degree() const;

			bool make_guess(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			bool propagate_na_constraints(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			virtual bool propagate_nb_constraints(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;

			bool complete_representative(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			bool minimally_extend_representative(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;

			virtual bool is_linearly_smaller(const function_in_out_t& minS) const;
			virtual bool is_linearly_minimal() const;
			virtual bool is_affinely_smaller(const function_in_out_t& minS) const;
			bool is_affinely_minimal() const;

			function_in_out_t get_linear_representative() const;
			function_in_out_t get_affine_representative() const;

			bool make_el_guess(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				linear_function_in_out_t& C,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				linear_function_in_out_t* minC = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			bool propagate_na_el_constraints(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				linear_function_in_out_t& C,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				linear_function_in_out_t* minC = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			bool propagate_nb_el_constraints(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				linear_function_in_out_t& C,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				linear_function_in_out_t* minC = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;

			bool complete_el_representative(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				linear_function_in_out_t& C,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				linear_function_in_out_t* minC = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;
			
			bool minimally_extend_el_representative(const function_in_out_t& minS,
				linear_permutation_in_t& A, linear_permutation_out_t& B,
				linear_function_in_out_t& C,
				bitset_in_t& nA, bitset_out_t& nB,
				function_in_out_t& R,
				linear_permutation_in_t* minA = NULL, linear_permutation_out_t* minB = NULL,
				linear_function_in_out_t* minC = NULL,
				function_in_out_t* minR = NULL,
				bool early_abort = true) const;

			bool is_extended_linearly_smaller(const function_in_out_t& minS) const;
			bool is_extended_linearly_minimal() const;
			bool is_extended_affinely_smaller(const function_in_out_t& minS) const;
			bool is_extended_affinely_minimal() const;

			void print_ddt() const;
			void print_bbt() const;

			std::string to_string(const char* preamble = "", const char* postamble = "", void (*sprintw)(char*, word_out_t) = &print_word_out, const size_t width = constants_out_t::WORD_CHARS, const char* wpreamble = "", const char* wpostamble = ",", const char* unknown = "*") const;
			std::string to_sage_string() const;
			std::string to_array_string() const;
	};

	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator<(const function_in_out_t& lhs, const function_in_out_t& rhs);
	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator>(const function_in_out_t& lhs, const function_in_out_t& rhs);
	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator==(const function_in_out_t& lhs, const function_in_out_t& rhs);
	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator!=(const function_in_out_t& lhs, const function_in_out_t& rhs);
	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator<=(const function_in_out_t& lhs, const function_in_out_t& rhs);
	template <size_t BITS_IN, size_t BITS_OUT>
	bool operator>=(const function_in_out_t& lhs, const function_in_out_t& rhs);

	template <size_t BITS_IN, size_t BITS_OUT>
	std::ostream& operator<<(std::ostream& os, const function_in_out_t& S);
}

using namespace std;
using namespace apn;

template <size_t BITS_IN, size_t BITS_OUT>
function_in_out_t::function(initializer_list<word_out_t> list) {
	const word_out_t* it = std::begin(list);
	for (word_in_t x = 0; x < list.size(); ++x) {
		set(x, *it);
		++it;
		if (x == constants_in_t::WORD_MAX)
			break;
	}
}

template <size_t BITS_IN, size_t BITS_OUT>
function_in_out_t::function(const array_in_out_t& a) {
	for (word_in_t x = 0; ; ++x) {
		set(x, a[x]);
		if (x == constants_in_t::WORD_MAX)
			break;
	}
}
template <size_t BITS_IN, size_t BITS_OUT>
std::string function_in_out_t::to_string(const char* preamble, const char* postamble, void (*sprintw)(char*, word_out_t), const size_t width, const char* wpreamble, const char* wpostamble, const char* unknown) const {
	const function_in_out_t &S = *this;
	const size_t wwidth = strlen(wpreamble) + width + strlen(wpostamble);
	const size_t bwidth = strlen(preamble) + wwidth*(1 << BITS_IN) + strlen(postamble) + 1; // +1 for null terminating character
	const size_t uwidth = strlen(unknown) + strlen(wpostamble);
	char* buffer = (char*) malloc(bwidth);
	string s;
	sprintf(buffer, "%s", preamble);
	for (word_in_t x = 0; ; ++x) {
		if (!(S.test(x))) {
			size_t i;
			for (i = 0; i < wwidth - uwidth; i++)
				sprintf(buffer + strlen(preamble) + wwidth*x + i, " ");
			sprintf(buffer + strlen(preamble) + wwidth*x + i, "%s%s", unknown, wpostamble);
		} else {
			sprintf(buffer + strlen(preamble) + wwidth*x, "%s", wpreamble);
			sprintw(buffer + strlen(preamble) + wwidth*x + strlen(wpreamble), S[x]);
			sprintf(buffer + strlen(preamble) + wwidth*x +  wwidth - strlen(wpostamble), "%s", wpostamble);
		}
		if (x == constants_in_t::WORD_MAX)
			break;
	}
	sprintf(buffer + strlen(preamble) + wwidth*(1 << BITS_IN), "%s", postamble);
	s = string(buffer);
	free(buffer);
	return s;
}

template <size_t BITS_IN, size_t BITS_OUT>
std::string function_in_out_t::to_sage_string() const {
	const char* preamble = "mq.SBox([";
	const char* postamble = "])";
	return to_string(preamble, postamble);
}

template <size_t BITS_IN, size_t BITS_OUT>
std::string function_in_out_t::to_array_string() const {
	const char* preamble = "{";
	const char* postamble = "}";
	return to_string(preamble, postamble, &print_hword_out, constants_out_t::HWORD_CHARS);
}

template <size_t BITS_IN, size_t BITS_OUT>
ostream& apn::operator<<(ostream& os, const function_in_out_t& S) {
	return os << S.to_string();
}	

/*
 * This methods only makes sense when the function_in_out_t is fully determined
 */
template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::print_ddt() const {
	const function_in_out_t& S = *this;
	const ddt_t& ddt = S.ddt;

	for (word_in_t x = 0; ; ++x) {
		for (word_out_t y = 0; ; ++y) {
			if (ddt[x][y])
				printf("%2d ", ddt[x][y]);
			else
				printf(" . ");

			if (y == constants_out_t::WORD_MAX)
				break;
		}
		printf("\n");

		if (x == constants_in_t::WORD_MAX)
			break;
	}
	printf("\n");
}

/*
 * This methods only makes sense when the function_in_out_t is fully determined
 */
template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::print_bbt() const {
	const function_in_out_t& S = *this;
	const bbt_t& bbt = S.bbt;

	for (word_in_t x = 0; ; ++x) {
		for (word_out_t y = 0; ; ++y) {
			if (bbt[x].test(y))
				printf("x ");
			else
				printf(". ");

			if (y == constants_out_t::WORD_MAX)
				break;
		}
		printf("\n");
		
		if (x == constants_in_t::WORD_MAX)
			break;
	}
	printf("\n");
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator<(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] < rhs[x])
				return true;
			if (lhs[x] > rhs[x])
				return false;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator>(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] > rhs[x])
				return true;
			if (lhs[x] < rhs[x])
				return false;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator<=(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] < rhs[x])
				return true;
			if (lhs[x] > rhs[x])
				return false;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return true;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator>=(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] > rhs[x])
				return true;
			if (lhs[x] < rhs[x])
				return false;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return true;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator==(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] != rhs[x])
				return false;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return true;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool apn::operator!=(const function_in_out_t& lhs, const function_in_out_t& rhs) {
	for (word_in_t x = 0; ; ++x) {
		if (lhs.test(x) && rhs.test(x)) {
			if (lhs[x] != rhs[x])
				return true;
			if (x == constants_in_t::WORD_MAX)
				break;
		} else {
			return false;
		}
	}
	return false;
}

/* Sets image of x to y and sets corresponding bits in [co]domain. */
template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set(const word_in_t x, const word_out_t y) {
	(*this)[x] = y;
	domain.set(x);
	codomain.set(y);
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::reset() {
	domain.reset();
	codomain.reset();
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::reset(word_in_t x) {
	if (test(x)) {
		word_out_t y = (*this)[x];
		domain.reset(x);
		codomain.reset(y);
		fix_codomain(y);
	}
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::fix_codomain(word_out_t y) {
	for (word_in_t x = 0; ; ++x) {
		if (test(x) && (*this)[x] == y) {
			codomain.set(y);
			return true;
		}
		if (x == constants_in_t::WORD_MAX)
			break;
	}
	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::test(word_in_t x) const {
	return domain.test(x);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::test_domain(word_in_t x) const {
	return test(x);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::test_codomain(word_out_t y) const {
	return codomain.test(y);
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set_domain(const bitset_in_t set)  {
	domain = set;
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set_codomain(const bitset_out_t set) {
	codomain = set;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_in_t function_in_out_t::get_domain() const {
	return domain;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_out_t function_in_out_t::get_codomain() const {
	return codomain;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_out_t function_in_out_t::image(bitset_in_t set) const {
	bitset_in_t range = get_domain() & set;
	bitset_out_t image;

	for (word_in_t x=0; ; ++x) {
		if (range[x])
			image.set((*this)[x]);
		if (x == constants_in_t::WORD_MAX)
			break;
	}

	return image;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_in_t function_in_out_t::preimages(word_out_t y) const {
	bitset_in_t preimages;

	if (test_codomain(y)) {
		for (word_in_t x = 0; ; ++x) {
			if (test(x) && (*this)[x] == y)
				preimages.set(x);
			if (x == constants_in_t::WORD_MAX)
				break;
		}
	}

	return preimages;
}

template <size_t BITS_IN, size_t BITS_OUT>
vector<word_in_t> function_in_out_t::preimages_vector(word_out_t y) const {
	vector<word_in_t> vec;

	for (word_in_t x = 0; ; ++x) {
		if (test(x) && (*this)[x] == y)
			vec.push_back(x);
		if (x == constants_in_t::WORD_MAX)
			break;
	}

	return vec;
}

template <size_t BITS_IN, size_t BITS_OUT>
bitset_in_t function_in_out_t::inverse_image(const bitset_out_t& set) const {
	bitset_out_t range = get_codomain() & set;
	bitset_in_t image;

	for (word_out_t y = 0; ; ++y) {
		if (range[y]) {
			image |= preimages(y);
		}
		if (y == constants_out_t::WORD_MAX)
			break;
	}/*x*/

	return image;
}

template <size_t BITS_IN, size_t BITS_OUT>
word_in_t function_in_out_t::get_position() const {
	return position;
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::inc_position() {
	position++;
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set_position(word_in_t pos){
	position = pos;
}

/*
** Check whether the partially left-defined Sbox S could be APN.
*/
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::extend_bbt(bool early_abort) {
	const function_in_out_t& S = *this;

	const word_out_t y = S[position];

	word_in_t diff_input;
	word_out_t diff_output;

	for (word_in_t x = 0; x < position; ++x) {

		diff_input  = x^position;
		diff_output = S[x]^y;

		/* If the BBT is already set here,
		 * that partial guessed permutation cannot be APN
		 */
		if (bbt[diff_input].test(diff_output)) {
			if (early_abort)
				return false;
			apn = false;
		}

		/* If not, update the BBT */
		bbt[diff_input].set(diff_output);

	}/*x*/

	return apn;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_apn() const {
	return apn;
}

template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::get_max_uniformity() const {
	return max_uniformity;
}

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set_max_uniformity(int uniformity) {
	max_uniformity = uniformity;
}

/*
 * This methods only makes sense when the function_in_out_t is fully determined
 */
template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::compute_ddt() {
	function_in_out_t& S = *this;

	int val;

	ddt = {{}}; // clean up previous computations
	for (word_in_t x = 0; ; ++x) {
		for (word_in_t diff_input = 0; ; ++diff_input) {
			val = ++ddt[diff_input][S[x]^S[x^diff_input]];
			if (diff_input != 0 && val > uniformity) {
				uniformity = val;
			}
			if (diff_input == constants_in_t::WORD_MAX)
				break;
		}/*a*/
		if (x == constants_in_t::WORD_MAX)
			break;
	}/*x*/

	return uniformity;
}

template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::get_uniformity() const {
	return uniformity;
}

/*
** Check whether the partially left-defined Sbox S could be APN.
*/
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::extend_ddt(bool early_abort) {
	const function_in_out_t& S = *this;

	const word_out_t y = S[position];

	word_in_t diff_input;
	word_out_t diff_output;

	for (word_in_t x = 0; x < position; ++x) {

		diff_input  = x^position;
		diff_output = S[x]^y;

		/* If the DDT already contains the max value here,
		 * that partial guessed permutation cannot be max-uniform.
		 */
		if (ddt[diff_input][diff_output] == max_uniformity) {
			if (early_abort) {
				return false;
			}
		}

		/* If not, update the DDT */
		ddt[diff_input][diff_output] += 2;

		/*
		if (x == constants_in_t::WORD_MAX)
			break;
		*/
	}/*x*/

	return true;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_involution() const {
	const function_in_out_t& S = *this;

	for (word_in_t x = 0; ; ++x) {
		if (!S.test(x)) {
			if (x == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		if (!S.test(S[x])) {
			if (x == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		if (x != S[S[x]])
			return false;
		if (x == constants_in_t::WORD_MAX)
			break;
	}/*x*/

	return true;
}

/******************* Affine Equivalence **************************************/
/* Return value is for early abortion */
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::make_guess(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	linear_permutation_in_t bak_A = A;
	linear_permutation_out_t bak_B = B;
	function_in_out_t bak_R = R;

	size_t dim;
	word_in_t x;

	/* Test all possible values, zero is not a valid value and will break stuff */
	for(word_in_t guess=1; ; ++guess) {

		/* Nothing to say if S is not defined there */
		if (!S.test(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* A is not a permutation_t */
		if (A.test_codomain(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* Pick the smallest unassigned element for A */
		dim = A.get_domain_dimension();
		x = 1 << dim;
#if DEBUG
		cout << "--------------------------------------\n";
		cout << "---- extending A ---------------------\n";
		cout << "--------------------------------------\n";
		cout << "S = " << S << endl;
		cout << "minS = " << minS << endl;
		cout << "A = " << A << endl;
		cout << "B = " << B << endl;
		cout << "R = " << R << endl;
		cout << "x = " << (long) x << "\nguess = " << (long) guess << "\n";
#endif
		/* Make the guess, no need to properly store A inverse */
		nA = A.quickly_add_refining_linear_span(x, guess, false);

		/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
		if (minimally_extend_representative(minS, A, B, nA, nB, R, minA, minB, minR, early_abort))
			return true;

		if (guess == constants_in_t::WORD_MAX)
			break;

		/* Then try next guess */
		A = bak_A;
		B = bak_B;
		R = bak_R;
		//nA = 0; // nA is set above
		nB = 0;

	} /* Loop on guesses */

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::propagate_na_constraints(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	word_in_t x;
	word_in_t y;
	word_out_t z;
	size_t dim;

#if SKIP_CLEANUP
#else
	bitset_out_t SAnA = S.image(A.image(nA));
	bitset_in_t common;
#endif

	/* A part */
	while (nA.any()) {
		x = ffs_in(nA);
		y = A[x];

		/* If S is not defined at A[x], nothing to say */
		if (!S.test(y)) {
			nA.reset(x);
			continue;
		}
	
		z = S[y];

		if (B.test_codomain(z)) {
			y = B.invert(z);
			R.set(x, y);
			nA.reset(x);
		} else {
			/* Next undefined value for B */
			dim = B.get_domain_dimension();
			y = 1 << dim;
			R.set(x, y);

#if SKIP_CLEANUP
			nB |= B.quickly_add_refining_linear_span(y, z, true);
#else
			nB |= B.get_domain() << y;
			common = B.add_refining_linear_span(y, z) & SAnA;
			/* S was maybe somehow linear */
			nB &= ~(B.inverse_image(common));
			nA &= ~(A.inverse_image(S.inverse_image(common)));
#endif
		}

		/* Early abort ?*/
		if (early_abort && R < minS) {
			return true;
		}
		/* No need to recurse if the result is larger than the current minimum */
		if (R > minS || (minR != NULL && R > *minR)) {
			return false;
		}
	}

	/* Now nA is empty and nB surely is not. */
	return minimally_extend_representative(minS, A, B, nA, nB, R, minA, minB, minR, early_abort);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::propagate_nb_constraints(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	linear_permutation_in_t bak_A;
	linear_permutation_out_t bak_B;
	function_in_out_t bak_R;
	bitset_out_t bak_nB;

	word_in_t x;
	word_out_t y;
	word_out_t z;
	size_t dim;

#if SKIP_CLEANUP
#else
	bitset_t common;
#endif

	/* B part */
	while (nB.any()) {
		y = ffs_out(nB);
		z = B[y];

#if SKIP_CLEANUP
		if (S.test_codomain(z) && (A.get_codomain() & S.preimages(z)).none())
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
	x = 1 << dim;

	/*
	 * Multiple choices for the preimage of z if S is not bijective.
	 * Worse, we have to check for smaller elements in the images of the vector spaces reaching z through S.
	 */
	word_out_t min_z = z;
	word_in_t u;
	word_out_t v;
	bitset_in_t preimages = S.preimages(z);
	bitset_in_t guesses = 0;
	for (word_in_t test = 0; ; test++) {
		if (!preimages.test(test)) {
			if (test == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		/* A[x] = test */
		for (word_in_t e = 0; e < x; e++) {
			u = test^A[e];
			if (!S.test(u))
				continue;
			v = S[u];
			if (!B.test_codomain(v))
				continue;
			if (v < min_z) {
				min_z = v;
				guesses = 0;
				guesses.set(u);
				continue;
			}
			if (v == min_z) {
				guesses.set(u);
				if (v == z)
					/* Avoid testing a preimage twice */
					preimages.reset(u);
			}
		}
		if (test == constants_in_t::WORD_MAX)
			break;
	}

	if (min_z != z)
		y = B.invert(min_z);
	R.set(x, y);

	/* Early abort ?*/
	if (early_abort && R < minS) {
		return true;
	}
	/* No need to recurse if the result is larger than the current minimum */
	if (R > minS || (minR != NULL && R > *minR)) {
		return false;
	}

	bak_A = A;
	bak_B = B;
	bak_nB = nB;
	bak_R = R;

	/*
	 * Now we have all preimages of the smallest element min_z
	 * in the vector spaces containing the preimages of z
	 */
	for (word_in_t guess = 0; ; guess++) {
		if (!guesses.test(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
#if DEBUG
		cout << "--------------------------------------\n";
		cout << "---- extending A from B---------------\n";
		cout << "--------------------------------------\n";
		cout << "S = " << S << endl;
		cout << "minS = " << minS << endl;
		cout << "A = " << A << endl;
		cout << "B = " << B << endl;
		cout << "R = " << R << endl;
		cout << "nA = " << nA << endl;
		cout << "nB = " << nA << endl;
		cout << "x = " << (long) x << endl;
		cout << "y = " << (long) y << endl;
		cout << "z = " << (long) z << endl;
		cout << "min_z = " << (long) min_z << endl;
		cout << "guess = " << (long) guess << endl;
#endif
		//cout << (long) x << " " << (long) guess << " " << (long) z << " " << (long) y << endl;

#if SKIP_CLEANUP
		nA = A.quickly_add_refining_linear_span(x, guess, false);
#else
		nA = A.get_domain() << x;
		common = S.image(A.add_refining_linear_span(x, guess));
		/* S was maybe somehow linear and we have to reset y anyway */
		nB &= ~(B.inverse_image(common));
		/* this could remove points in nA where S is not defined */
		//nA &= ~(A.inverse_image(S.inverse_image(common)));
#endif

		/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
		if (minimally_extend_representative(minS, A, B, nA, nB, R, minA, minB, minR, early_abort))
		//if (propagate_na_constraints(minS, A, B, nA, nB, R, minA, minB, minR, early_abort))
			return true;

		if (guess == constants_in_t::WORD_MAX)
			break;

		/* Then try next guess */
		A = bak_A;
		B = bak_B;
		R = bak_R;
		//nA = 0;
		nB = bak_nB;
	} /* Try next possible inverse */

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::complete_representative(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	/* Complete R */
	size_t dim = A.get_domain_dimension();
	//linear_permutation_t Binv = B.dirty_inverse();
	for (word_in_t x = 0; x < (1 << dim); ++x) {
		word_in_t y = A[x];
		if (!R.test(x) && S.test(y)) {
			word_out_t z = S[y];
			if (B.test_codomain(z)) {
				R.set(x, B.invert(z));
			}
			/*if (Binv.test(z)) {
				R.set(x, Binv[z]);
			}*/
		}
		if (x == constants_in_t::WORD_MAX)
			break;
	}

	if (early_abort && R < minS)
		return true;

	/* Did we find a smaller representative? */
	if (minR != NULL && !(*minR < R)) {
		*minR = R;
		*minA = A;
		*minB = B;
	}

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::minimally_extend_representative(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	/* The caller can have set something in nA, treat it */
	if (nA.any()) {
		return propagate_na_constraints(minS, A, B, nA, nB, R, minA, minB, minR, early_abort);
	}

	/* nA was empty but nB is not. */
	if (nB.any()) {
		return propagate_nb_constraints(minS, A, B, nA, nB, R, minA, minB, minR, early_abort);
	}

	/* Are finished? */
	if (((~A.get_codomain()) & S.get_domain()).none()) {
		//if (nA.any() || nB.any()) abort(); // something wrong happened !?!
		return complete_representative(minS, A, B, R, minA, minB, minR, early_abort);
	}

	return make_guess(minS, A, B, nA, nB, R, minA, minB, minR, early_abort);
}

/* Assert S[0] is set */
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_linearly_smaller(const function_in_out_t& minS) const {
	const function_in_out_t& S = *this;
	linear_permutation_in_t A;
	linear_permutation_out_t B;
	function_in_out_t R;

	if (!test(0) || !minS.test(0))
		return false;

	//A.set(0, 0);
	//B.set(0, 0);
	bitset_in_t nA = 0;
	bitset_out_t nB = 0;
	if (S[0] != 0) {
		if (minS[0] == 0)
			return false;
		if (minS[0] > 1)
			return true;
		nA.set(0);
		nB.set(0); 
	} else {
		/* S[0] == 0 */
		if (minS[0] != 0)
			return true;
		if (minS.test(1) && minS[1] > 1)
			return true;
		R.set(0, 0);
	}

	return minimally_extend_representative(minS, A, B, nA, nB, R, NULL, NULL, NULL, true);
}


template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_linearly_minimal() const {
	const function_in_out_t& S = *this;

	if (test(0)) {
		if (S[0] != 0) {
			if (S[0] != 1)
				return false;
		} else {
			if (test(1) && S[1] > 1)
				return false;
		}
	}

	return !is_linearly_smaller(S);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_affinely_smaller(const function_in_out_t& minS) const {
	const function_in_out_t& S = *this;

	if (minS.test(0) && minS[0] != 0)
		return true;

	/* Iterate through shifts of the original S-box with 0 in zero position */
	for (word_in_t a = 0; ; ++a) {
		if (!test(a)) {
			if (a == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		function_in_out_t shifted_function_in_out_t;
		/* Set shifted[0] = 0 */
		word_out_t b = S[a];

		for (word_in_t x = 0; ; ++x) {
			word_out_t y = x^a;
			if (test(y)) {
				shifted_function_in_out_t.set(x, S[y]^b);
			}
			if (x == constants_in_t::WORD_MAX)
				break;
		}/*x*/

		if (shifted_function_in_out_t.is_linearly_smaller(minS)) {
			return true;
		}

		if (a == constants_in_t::WORD_MAX)
			break;
	}/*a*/

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_affinely_minimal() const {
	const function_in_out_t& S = *this;

	if (test(0) && S[0] != 0) {
		return false;
	}

	if (test(1) && S[1] > 1) {
		return false;
	}

	return !is_affinely_smaller(S);
}

/* Assert S[0] is set */
template <size_t BITS_IN, size_t BITS_OUT>
function_in_out_t function_in_out_t::get_linear_representative() const {

	const function_in_out_t& S = *this;

	linear_permutation_in_t A;
	linear_permutation_out_t B;
	function_in_out_t R;

	linear_permutation_in_t minA;
	linear_permutation_out_t minB;
	function_in_out_t minR;

	bitset_in_t nA = 0;
	bitset_out_t nB = 0;

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
template <size_t BITS_IN, size_t BITS_OUT>
function_in_out_t function_in_out_t::get_affine_representative() const {

	const function_in_out_t& S = *this;
	function_in_out_t minS = S;

	for (word_in_t a = 0; a<=constants_in_t::WORD_MAX; ++a) {
		for(word_out_t b = 0; b<=constants_out_t::WORD_MAX; ++b) {

			function_in_out_t shifted_function_in_out_t;

			for (word_in_t x = 0; x<=constants_in_t::WORD_MAX; ++x) {
				shifted_function_in_out_t.set(x, S[x^a]^b);
			}/*x*/
			
			function_in_out_t tmp = shifted_function_in_out_t.get_linear_representative();
			if(tmp<minS) minS = tmp;
			
		}/*b*/
	}/*a*/

	return minS;
}

/******************* Extended Affine Equivalence *****************************/
/* Just the same as above... */
/* Return value is for early abortion */
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::make_el_guess(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	linear_function_in_out_t& C,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	linear_function_in_out_t* minC,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	linear_permutation_in_t bak_A = A;
	linear_permutation_out_t bak_B = B;
	linear_function_in_out_t bak_C = C;
	function_in_out_t bak_R = R;

	size_t dim;
	word_in_t x;

	/* Test all possible values, zero is not a valid value and will break stuff */
	for(word_in_t guess = 1; ; ++guess) {

		/* Nothing to say if S is not defined there */
		if (!S.test(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* A is not a permutation_t */
		if (A.test_codomain(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}

		/* Pick the smallest unassigned element for A */
		dim = A.get_domain_dimension();
		x = 1 << dim;
#if DEBUG
		cout << "--------------------------------------\n";
		cout << "---- extending A ---------------------\n";
		cout << "--------------------------------------\n";
		cout << "S = " << S << endl;
		cout << "minS = " << minS << endl;
		cout << "A = " << A << endl;
		cout << "B = " << B << endl;
		cout << "C = " << C << endl;
		cout << "R = " << R << endl;
		cout << "x = " << (long) x << "\nguess = " << (long) guess << "\n";
#endif
		/* Make the guess */
		nA = A.quickly_add_refining_linear_span(x, guess, false);


		/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
		if (minimally_extend_el_representative(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort))
			return true;

		if (guess == constants_in_t::WORD_MAX)
			break;

		/* Then try next guess */
		A = bak_A;
		B = bak_B;
		C = bak_C;
		R = bak_R;
		//nA = 0; // nA is set above
		nB = 0;

	} /* Loop on guesses */

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::propagate_na_el_constraints(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	linear_function_in_out_t& C,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	linear_function_in_out_t* minC,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	word_in_t x;
	word_in_t y;
	word_out_t z;

#if SKIP_CLEANUP
#else
	bitset_out_t SAnA = S.image(A.image(nA));
	bitset_in_t common;
#endif

	/* A part */
	while (nA.any()) {
		x = ffs_in(nA);
		y = A[x];

		/* If S is not defined at A[x], nothing to say */
		if (!S.test(y)) {
			nA.reset(x);
			continue;
		}
	
		z = S[y];

		if (B.test_codomain(z)) {
			y = B.invert(z);
			if (C.test(x)) {
				R.set(x, y^C[x]);
			} else if (__builtin_popcount(x) == 1) {
				C.add_linear_span(x, y);
				R.set(x, 0);
			}
			nA.reset(x);
		} else {
			break;
		}

		/* Early abort ?*/
		if (early_abort && R < minS) {
			return true;
		}
		/* No need to recurse if the result is larger than the current minimum */
		if (R > minS || (minR != NULL && R > *minR)) {
			return false;
		}
	}

	if (nA.any()) {
		linear_permutation_in_t bak_A = A;
		linear_permutation_out_t bak_B = B;
		linear_function_in_out_t bak_C = C;
		bitset_in_t bak_nA = nA;
		bitset_out_t bak_nB = nB;
		function_in_out_t bak_R = R;

		/* Make a guess for B and ensure EA minimality through C */
		for (word_out_t guess = 1; ; ++guess) {

			/* Get undefined value for B */
			if (B.test(guess)) {
				if (guess == constants_out_t::WORD_MAX)
					break;
				else
					continue;
			}

			if (__builtin_popcount(x) == 1) {
				C.add_linear_span(x, guess);
				R.set(x, 0);
			}

#if SKIP_CLEANUP
			nB |= B.add_linear_span(guess, z);
#else
#error "To the infinity and beyond."
#endif

			if (minimally_extend_el_representative(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort))
				return true;

			if (guess == constants_out_t::WORD_MAX)
				break;

			/* Then try next guess */
			A = bak_A;
			B = bak_B;
			C = bak_C;
			R = bak_R;
			nA = bak_nA;
			nB = bak_nB;
		} /* Loop on guesses */

		/* Recursion found nothing smaller */
		return false;
	}

	/* No guess where made and nA is empty. */
	return minimally_extend_el_representative(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::propagate_nb_el_constraints(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	linear_function_in_out_t& C,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	linear_function_in_out_t* minC,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	linear_permutation_in_t bak_A;
	linear_permutation_out_t bak_B;
	linear_function_in_out_t bak_C;
	function_in_out_t bak_R;
	bitset_out_t bak_nB;

	word_in_t x;
	word_out_t y;
	word_out_t z;
	size_t dim;

#if SKIP_CLEANUP
#else
	bitset_t common;
#endif

	/* B part */
	while (nB.any()) {
		y = ffs_out(nB);
		z = B[y];

#if SKIP_CLEANUP
		if (S.test_codomain(z) && (A.get_codomain() & S.preimages(z)).none())
#else
		if (S.test_codomain(z))
#endif
		{
			break;
		}
		nB.reset(y);
	}

	if (nB.none())
		return minimally_extend_el_representative(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort);

	/* Next undefined value for A*/
	dim = A.get_domain_dimension();
	x = 1 << dim;

	/*
	 * Multiple choices for the preimage of z if S is not bijective.
	 * Worse, we have to check for elements in the images of the vector spaces reaching z through S.
	 */
	word_in_t u;
	word_out_t v;
	bitset_in_t preimages = S.preimages(z);
	bitset_in_t guesses = 0;
	for (word_in_t test = 0; ; test++) {
		if (!preimages.test(test)) {
			if (test == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		/* A[x] = test */
		for (word_in_t e = 0; e < x; e++) {
			u = test^A[e];
			if (!S.test(u))
				continue;
			v = S[u];
			if (!B.test_codomain(v))
				continue;
			guesses.set(u);
			preimages.reset(u);
		}
		if (test == constants_in_t::WORD_MAX)
			break;
	}


	bak_A = A;
	bak_B = B;
	bak_C = C;
	bak_nB = nB;
	bak_R = R;

	/*
	 * Now we have all preimages of the elements
	 * in the vector spaces containing the preimages of z
	 */
	for (word_in_t guess = 0; ; guess++) {
		if (!guesses.test(guess)) {
			if (guess == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
#if DEBUG
		cout << "--------------------------------------\n";
		cout << "---- extending A from B---------------\n";
		cout << "--------------------------------------\n";
		cout << "S = " << S << endl;
		cout << "minS = " << minS << endl;
		cout << "A = " << A << endl;
		cout << "B = " << B << endl;
		cout << "C = " << C << endl;
		cout << "R = " << R << endl;
		cout << "nA = " << nA << endl;
		cout << "nB = " << nA << endl;
		cout << "x = " << (long) x << endl;
		cout << "y = " << (long) y << endl;
		cout << "z = " << (long) z << endl;
		cout << "min_z = " << (long) min_z << endl;
		cout << "guess = " << (long) guess << endl;
#endif

		y = B.invert(S[guess]);
		if (C.test(x)) {
			R.set(x, y^C[x]);
		} else if (__builtin_popcount(x) == 1) {
			C.add_linear_span(x, y);
			R.set(x, 0);
		}

		/* Early abort ?*/
		if (early_abort && R < minS) {
			return true;
		}
		/* No need to recurse if the result is larger than the current minimum */
		if (R > minS || (minR != NULL && R > *minR)) {
			return false;
		}

#if SKIP_CLEANUP
		nA = A.quickly_add_refining_linear_span(x, guess, false);
		//nA = A.quickly_add_linear_span(x, guess);
#else
		nA = A.get_domain() << x;
		common = S.image(A.add_refining_linear_span(x, guess));
		/* S was maybe somehow linear and we have to reset y anyway */
		nB &= ~(B.inverse_image(common));
		/* this could remove points in nA where S is not defined */
		//nA &= ~(A.inverse_image(S.inverse_image(common)));
#endif

		/* nA is non empty so this will call propagate_na_constraints anyway but make recursion more intricate */
		if (minimally_extend_el_representative(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort))
			return true;

		if (guess == constants_in_t::WORD_MAX)
			break;

		/* Then try next guess */
		A = bak_A;
		B = bak_B;
		C = bak_C;
		R = bak_R;
		//nA = 0;
		nB = bak_nB;
	} /* Try next possible inverse */

	return false;
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::complete_el_representative(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	linear_function_in_out_t& C,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	linear_function_in_out_t* minC,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	/* Complete R */
	size_t dim = A.get_domain_dimension();
	//linear_permutation_t Binv = B.dirty_inverse();
	for (word_in_t x = 0; x < (1 << dim); ++x) {
		word_in_t y = A[x];
		if (!R.test(x) && S.test(y)) {
			word_out_t z = S[y];
			if (B.test_codomain(z)) {
				y = B.invert(z);
				if (C.test(x)) {
					R.set(x, y^C[x]);
				} else if (__builtin_popcount(x) == 1) {
					C.add_linear_span(x, y);
					R.set(x, 0);
				}
			}
			/*if (Binv.test(z)) {
				R.set(x, Binv[z]);
			}*/
		}
		if (x == constants_in_t::WORD_MAX)
			break;
	}

	if (early_abort && R < minS)
		return true;

	/* Did we find a smaller representative? */
	if (minR != NULL && !(*minR < R)) {
		*minR = R;
		*minA = A;
		*minB = B;
		*minC = C;
	}

	return false;
}


template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::minimally_extend_el_representative(
	const function_in_out_t& minS,
	linear_permutation_in_t& A, linear_permutation_out_t& B,
	linear_function_in_out_t& C,
	bitset_in_t& nA, bitset_out_t& nB,
	function_in_out_t& R,
	linear_permutation_in_t* minA, linear_permutation_out_t* minB,
	linear_function_in_out_t* minC,
	function_in_out_t* minR,
	bool early_abort
) const {
	const function_in_out_t& S = *this;

	/* The caller can have set something in nA, treat it */
	if (nA.any()) {
		return propagate_na_el_constraints(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort);
	}

	/* nA was empty but nB is not. */
	if (nB.any()) {
		return propagate_nb_el_constraints(minS, A, B, C, nA, nB, R, minA, minB,minC, minR, early_abort);
	}

	/* Are finished? */
	if (((~A.get_codomain()) & S.get_domain()).none()) {
		//if (nA.any() || nB.any()) abort(); // something wrong happened !?!
		return complete_el_representative(minS, A, B, C, R, minA, minB, minC, minR, early_abort);
	}

	return make_el_guess(minS, A, B, C, nA, nB, R, minA, minB, minC, minR, early_abort);
}

/* Assert S[0] is set */
template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_extended_linearly_smaller(const function_in_out_t& minS) const {
	const function_in_out_t& S = *this;
	linear_permutation_in_t A;
	linear_permutation_out_t B;
	linear_function_in_out_t C;
	function_in_out_t R;

	if (!test(0) || !minS.test(0))
		return false;

	//A.set(0, 0);
	//B.set(0, 0);
	bitset_in_t nA = 0;
	bitset_out_t nB = 0;
	if (S[0] != 0) {
		if (minS[0] == 0)
			return false;
		if (minS[0] > 1)
			return true;
		nA.set(0);
		nB.set(0); 
	} else {
		/* S[0] == 0 */
		if (minS[0] != 0)
			return true;
		if (minS.test(1) && minS[1] > 1)
			return true;
		R.set(0, 0);
	}

	return minimally_extend_el_representative(minS, A, B, C, nA, nB, R, NULL, NULL, NULL, NULL, true);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_extended_linearly_minimal() const {
	const function_in_out_t& S = *this;

	if (test(0)) {
		if (S[0] != 0) {
			if (S[0] != 1)
				return false;
		} else {
			if (test(1) && S[1] > 1)
				return false;
		}
	}

	return !is_extended_linearly_smaller(S);
}

template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_extended_affinely_smaller(const function_in_out_t& minS) const {
	const function_in_out_t& S = *this;

	if (minS.test(0) && minS[0] != 0)
		return true;

	/* Iterate through shifts of the original S-box with 0 in zero position */
	for (word_in_t a = 0; ; ++a) {
		if (!test(a)) {
			if (a == constants_in_t::WORD_MAX)
				break;
			else
				continue;
		}
		function_in_out_t shifted_function_in_out_t;
		/* Set shifted[0] = 0 */
		word_out_t b = S[a];

		for (word_in_t x = 0; ; ++x) {
			word_out_t y = x^a;
			if (test(y)) {
				shifted_function_in_out_t.set(x, S[y]^b);
			}
			if (x == constants_in_t::WORD_MAX)
				break;
		}/*x*/

		if (shifted_function_in_out_t.is_extended_linearly_smaller(minS)) {
			return true;
		}

		if (a == constants_in_t::WORD_MAX)
			break;
	}/*a*/

	return false;
}


template <size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::is_extended_affinely_minimal() const {
	const function_in_out_t& S = *this;

	if (test(0) && S[0] != 0) {
		return false;
	}

	for (word_in_t pos = 1; pos <= (1 << (BITS_IN-1)); pos <<= 1) {
		if (test(pos) && S[pos] != 0) {
			return false;
		}
	}

	return !is_extended_affinely_smaller(S);
}

/******************* Degree **************************************************/

template <size_t BITS_IN, size_t BITS_OUT>
void function_in_out_t::set_max_degree(int degree) {
	max_degree = degree;
}

template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::get_max_degree() const {
	return max_degree;
}

/*
 * Check whether the values in S[0..position-1] can be used to derive S[position]
 * when the function_in_out_t S is assumed to be of known degree
 */
template<size_t BITS_IN, size_t BITS_OUT>
bool function_in_out_t::propagate_degree_constraints(__attribute__((__unused__)) word_in_t position) {
	return false;
}

namespace apn{
#define BITS_IN 5
#define BITS_OUT 5
template <>
bool function_in_out_t::propagate_degree_constraints(word_in_t position) {
	function_in_out_t& S = *this;
	if (max_degree == 2) {
		/* 5-bit degree-2 */
		switch (position) { 
			case  7: S.set(position, S[0]^S[1]^S[2]^S[3]^S[4]^S[5]^S[6]); break;
			case 11: S.set(position, S[0]^S[1]^S[2]^S[3]^S[8]^S[9]^S[10]); break;
			case 13: S.set(position, S[0]^S[1]^S[4]^S[5]^S[8]^S[9]^S[12]); break;
			case 14: S.set(position, S[0]^S[2]^S[4]^S[6]^S[8]^S[10]^S[12]); break;
			case 15: S.set(position, S[0]^S[1]^S[2]^S[3]^S[12]^S[13]^S[14]); break;
			case 19: S.set(position, S[0]^S[1]^S[2]^S[3]^S[16]^S[17]^S[18]); break;
			case 21: S.set(position, S[0]^S[1]^S[4]^S[5]^S[16]^S[17]^S[20]); break;
			case 22: S.set(position, S[0]^S[2]^S[4]^S[6]^S[16]^S[18]^S[20]); break;
			case 23: S.set(position, S[0]^S[1]^S[2]^S[3]^S[20]^S[21]^S[22]); break;
			case 25: S.set(position, S[0]^S[1]^S[8]^S[9]^S[16]^S[17]^S[24]); break;
			case 26: S.set(position, S[0]^S[2]^S[8]^S[10]^S[16]^S[18]^S[24]); break;
			case 27: S.set(position, S[0]^S[1]^S[2]^S[3]^S[24]^S[25]^S[26]); break;
			case 28: S.set(position, S[0]^S[4]^S[8]^S[12]^S[16]^S[20]^S[24]); break;
			case 29: S.set(position, S[0]^S[1]^S[4]^S[5]^S[24]^S[25]^S[28]); break;
			case 30: S.set(position, S[0]^S[2]^S[4]^S[6]^S[24]^S[26]^S[28]); break;
			case 31: S.set(position, S[0]^S[1]^S[2]^S[3]^S[28]^S[29]^S[30]); break;
			default: return false;
		}
	} else if (max_degree == 3) {
		/* 5-bit degree-3 */
		switch (position) {
			case 15: S.set(position, S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 23: S.set(position, S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 27: S.set(position, S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 29: S.set(position, S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 30: S.set(position, S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 31: S.set(position, S[28]^S[26]^S[25]^S[22]^S[21]^S[19]^S[16]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1]); break;
			default: return false;
		}
	} else {
		return false;
	}

	return true;
}
#undef BITS_IN
#undef BITS_OUT

#define BITS_IN 6
#define BITS_OUT 4
template <>
bool function_in_out_t::propagate_degree_constraints(word_in_t position) {
	function_in_out_t& S = *this;

	if (max_degree == 2) {
		/* 6-bit degree-2 */
		/* Note: It's already proved that there is no degree-2 6-bit APN */
		switch (position) {
			case  7: S.set(position, S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 11: S.set(position, S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 13: S.set(position, S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 14: S.set(position, S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 15: S.set(position, S[12]^S[10]^S[9]^S[6]^S[5]^S[3]^S[0]); break;
			case 19: S.set(position, S[18]^S[17]^S[16]^S[3]^S[2]^S[1]^S[0]); break;
			case 21: S.set(position, S[20]^S[17]^S[16]^S[5]^S[4]^S[1]^S[0]); break;
			case 22: S.set(position, S[20]^S[18]^S[16]^S[6]^S[4]^S[2]^S[0]); break;
			case 23: S.set(position, S[20]^S[18]^S[17]^S[6]^S[5]^S[3]^S[0]); break;
			case 25: S.set(position, S[24]^S[17]^S[16]^S[9]^S[8]^S[1]^S[0]); break;
			case 26: S.set(position, S[24]^S[18]^S[16]^S[10]^S[8]^S[2]^S[0]); break;
			case 27: S.set(position, S[24]^S[18]^S[17]^S[10]^S[9]^S[3]^S[0]); break;
			case 28: S.set(position, S[24]^S[20]^S[16]^S[12]^S[8]^S[4]^S[0]); break;
			case 29: S.set(position, S[24]^S[20]^S[17]^S[12]^S[9]^S[5]^S[0]); break;
			case 30: S.set(position, S[24]^S[20]^S[18]^S[12]^S[10]^S[6]^S[0]); break;
			case 31: S.set(position, S[24]^S[20]^S[19]^S[12]^S[11]^S[7]^S[0]); break;
			case 35: S.set(position, S[34]^S[33]^S[32]^S[3]^S[2]^S[1]^S[0]); break;
			case 37: S.set(position, S[36]^S[33]^S[32]^S[5]^S[4]^S[1]^S[0]); break;
			case 38: S.set(position, S[36]^S[34]^S[32]^S[6]^S[4]^S[2]^S[0]); break;
			case 39: S.set(position, S[36]^S[34]^S[33]^S[6]^S[5]^S[3]^S[0]); break;
			case 41: S.set(position, S[40]^S[33]^S[32]^S[9]^S[8]^S[1]^S[0]); break;
			case 42: S.set(position, S[40]^S[34]^S[32]^S[10]^S[8]^S[2]^S[0]); break;
			case 43: S.set(position, S[40]^S[34]^S[33]^S[10]^S[9]^S[3]^S[0]); break;
			case 44: S.set(position, S[40]^S[36]^S[32]^S[12]^S[8]^S[4]^S[0]); break;
			case 45: S.set(position, S[40]^S[36]^S[33]^S[12]^S[9]^S[5]^S[0]); break;
			case 46: S.set(position, S[40]^S[36]^S[34]^S[12]^S[10]^S[6]^S[0]); break;
			case 47: S.set(position, S[40]^S[36]^S[35]^S[12]^S[11]^S[7]^S[0]); break;
			case 49: S.set(position, S[48]^S[33]^S[32]^S[17]^S[16]^S[1]^S[0]); break;
			case 50: S.set(position, S[48]^S[34]^S[32]^S[18]^S[16]^S[2]^S[0]); break;
			case 51: S.set(position, S[48]^S[34]^S[33]^S[18]^S[17]^S[3]^S[0]); break;
			case 52: S.set(position, S[48]^S[36]^S[32]^S[20]^S[16]^S[4]^S[0]); break;
			case 53: S.set(position, S[48]^S[36]^S[33]^S[20]^S[17]^S[5]^S[0]); break;
			case 54: S.set(position, S[48]^S[36]^S[34]^S[20]^S[18]^S[6]^S[0]); break;
			case 55: S.set(position, S[48]^S[36]^S[35]^S[20]^S[19]^S[7]^S[0]); break;
			case 56: S.set(position, S[48]^S[40]^S[32]^S[24]^S[16]^S[8]^S[0]); break;
			case 57: S.set(position, S[48]^S[40]^S[33]^S[24]^S[17]^S[9]^S[0]); break;
			case 58: S.set(position, S[48]^S[40]^S[34]^S[24]^S[18]^S[10]^S[0]); break;
			case 59: S.set(position, S[48]^S[40]^S[35]^S[24]^S[19]^S[11]^S[0]); break;
			case 60: S.set(position, S[48]^S[40]^S[36]^S[24]^S[20]^S[12]^S[0]); break;
			case 61: S.set(position, S[48]^S[40]^S[37]^S[24]^S[21]^S[13]^S[0]); break;
			case 62: S.set(position, S[48]^S[40]^S[38]^S[24]^S[22]^S[14]^S[0]); break;
			case 63: S.set(position, S[48]^S[40]^S[39]^S[24]^S[23]^S[15]^S[0]); break;
			default: return false;
		}
	} else if (max_degree == 3) {
		/* 6-bit degree-3 */
		switch (position) {
			case 15: S.set(position, S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 23: S.set(position, S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 27: S.set(position, S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 29: S.set(position, S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 30: S.set(position, S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 31: S.set(position, S[28]^S[26]^S[25]^S[22]^S[21]^S[19]^S[16]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1]); break;
			case 39: S.set(position, S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 43: S.set(position, S[42]^S[41]^S[40]^S[35]^S[34]^S[33]^S[32]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 45: S.set(position, S[44]^S[41]^S[40]^S[37]^S[36]^S[33]^S[32]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 46: S.set(position, S[44]^S[42]^S[40]^S[38]^S[36]^S[34]^S[32]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 47: S.set(position, S[44]^S[42]^S[41]^S[38]^S[37]^S[35]^S[32]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1]); break;
			case 51: S.set(position, S[50]^S[49]^S[48]^S[35]^S[34]^S[33]^S[32]^S[19]^S[18]^S[17]^S[16]^S[3]^S[2]^S[1]^S[0]); break;
			case 53: S.set(position, S[52]^S[49]^S[48]^S[37]^S[36]^S[33]^S[32]^S[21]^S[20]^S[17]^S[16]^S[5]^S[4]^S[1]^S[0]); break;
			case 54: S.set(position, S[52]^S[50]^S[48]^S[38]^S[36]^S[34]^S[32]^S[22]^S[20]^S[18]^S[16]^S[6]^S[4]^S[2]^S[0]); break;
			case 55: S.set(position, S[52]^S[50]^S[49]^S[38]^S[37]^S[35]^S[32]^S[22]^S[21]^S[19]^S[16]^S[7]^S[4]^S[2]^S[1]); break;
			case 57: S.set(position, S[56]^S[49]^S[48]^S[41]^S[40]^S[33]^S[32]^S[25]^S[24]^S[17]^S[16]^S[9]^S[8]^S[1]^S[0]); break;
			case 58: S.set(position, S[56]^S[50]^S[48]^S[42]^S[40]^S[34]^S[32]^S[26]^S[24]^S[18]^S[16]^S[10]^S[8]^S[2]^S[0]); break;
			case 59: S.set(position, S[56]^S[50]^S[49]^S[42]^S[41]^S[35]^S[32]^S[26]^S[25]^S[19]^S[16]^S[11]^S[8]^S[2]^S[1]); break;
			case 60: S.set(position, S[56]^S[52]^S[48]^S[44]^S[40]^S[36]^S[32]^S[28]^S[24]^S[20]^S[16]^S[12]^S[8]^S[4]^S[0]); break;
			case 61: S.set(position, S[56]^S[52]^S[49]^S[44]^S[41]^S[37]^S[32]^S[28]^S[25]^S[21]^S[16]^S[13]^S[8]^S[4]^S[1]); break;
			case 62: S.set(position, S[56]^S[52]^S[50]^S[44]^S[42]^S[38]^S[32]^S[28]^S[26]^S[22]^S[16]^S[14]^S[8]^S[4]^S[2]); break;
			case 63: S.set(position, S[56]^S[52]^S[51]^S[44]^S[43]^S[39]^S[32]^S[28]^S[27]^S[23]^S[16]^S[15]^S[8]^S[4]^S[3]); break;
			default: return false;
		}
	} else if (max_degree == 4) {
		/* 6-bit degree-4 */
		switch (position) {
			case 31: S.set(position, S[30]^S[29]^S[28]^S[27]^S[26]^S[25]^S[24]^S[23]^S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[15]^S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 47: S.set(position, S[46]^S[45]^S[44]^S[43]^S[42]^S[41]^S[40]^S[39]^S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[15]^S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 55: S.set(position, S[54]^S[53]^S[52]^S[51]^S[50]^S[49]^S[48]^S[39]^S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[23]^S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 59: S.set(position, S[58]^S[57]^S[56]^S[51]^S[50]^S[49]^S[48]^S[43]^S[42]^S[41]^S[40]^S[35]^S[34]^S[33]^S[32]^S[27]^S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 61: S.set(position, S[60]^S[57]^S[56]^S[53]^S[52]^S[49]^S[48]^S[45]^S[44]^S[41]^S[40]^S[37]^S[36]^S[33]^S[32]^S[29]^S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 62: S.set(position, S[60]^S[58]^S[56]^S[54]^S[52]^S[50]^S[48]^S[46]^S[44]^S[42]^S[40]^S[38]^S[36]^S[34]^S[32]^S[30]^S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 63: S.set(position, S[60]^S[58]^S[57]^S[54]^S[53]^S[51]^S[48]^S[46]^S[45]^S[43]^S[40]^S[39]^S[36]^S[34]^S[33]^S[30]^S[29]^S[27]^S[24]^S[23]^S[20]^S[18]^S[17]^S[15]^S[12]^S[10]^S[9]^S[6]^S[5]^S[3]^S[0]); break;
			default: return false;
		}
	} else {
		return false;
	}

	return true;
}
#undef BITS_IN
#undef BITS_OUT

#define BITS_IN 7
#define BITS_OUT 7
template <>
bool function_in_out_t::propagate_degree_constraints(word_in_t position) {
	function_in_out_t& S = *this;

	if (max_degree == 2) {
		/* 7-bit degree-2 */
		switch (position) {
			case  7: S.set(position, S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0]); break;
			case 11: S.set(position, S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0]); break;
			case 13: S.set(position, S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0]); break;
			case 14: S.set(position, S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0]); break;
			case 15: S.set(position, S[12]^S[10]^S[9]^S[6]^S[5]^S[3]^S[0]); break;
			case 19: S.set(position, S[18]^S[17]^S[16]^S[3]^S[2]^S[1]^S[0]); break;
			case 21: S.set(position, S[20]^S[17]^S[16]^S[5]^S[4]^S[1]^S[0]); break;
			case 22: S.set(position, S[20]^S[18]^S[16]^S[6]^S[4]^S[2]^S[0]); break;
			case 23: S.set(position, S[20]^S[18]^S[17]^S[6]^S[5]^S[3]^S[0]); break;
			case 25: S.set(position, S[24]^S[17]^S[16]^S[9]^S[8]^S[1]^S[0]); break;
			case 26: S.set(position, S[24]^S[18]^S[16]^S[10]^S[8]^S[2]^S[0]); break;
			case 27: S.set(position, S[24]^S[18]^S[17]^S[10]^S[9]^S[3]^S[0]); break;
			case 28: S.set(position, S[24]^S[20]^S[16]^S[12]^S[8]^S[4]^S[0]); break;
			case 29: S.set(position, S[24]^S[20]^S[17]^S[12]^S[9]^S[5]^S[0]); break;
			case 30: S.set(position, S[24]^S[20]^S[18]^S[12]^S[10]^S[6]^S[0]); break;
			case 31: S.set(position, S[24]^S[20]^S[19]^S[12]^S[11]^S[7]^S[0]); break;
			case 35: S.set(position, S[34]^S[33]^S[32]^S[3]^S[2]^S[1]^S[0]); break;
			case 37: S.set(position, S[36]^S[33]^S[32]^S[5]^S[4]^S[1]^S[0]); break;
			case 38: S.set(position, S[36]^S[34]^S[32]^S[6]^S[4]^S[2]^S[0]); break;
			case 39: S.set(position, S[36]^S[34]^S[33]^S[6]^S[5]^S[3]^S[0]); break;
			case 41: S.set(position, S[40]^S[33]^S[32]^S[9]^S[8]^S[1]^S[0]); break;
			case 42: S.set(position, S[40]^S[34]^S[32]^S[10]^S[8]^S[2]^S[0]); break;
			case 43: S.set(position, S[40]^S[34]^S[33]^S[10]^S[9]^S[3]^S[0]); break;
			case 44: S.set(position, S[40]^S[36]^S[32]^S[12]^S[8]^S[4]^S[0]); break;
			case 45: S.set(position, S[40]^S[36]^S[33]^S[12]^S[9]^S[5]^S[0]); break;
			case 46: S.set(position, S[40]^S[36]^S[34]^S[12]^S[10]^S[6]^S[0]); break;
			case 47: S.set(position, S[40]^S[36]^S[35]^S[12]^S[11]^S[7]^S[0]); break;
			case 49: S.set(position, S[48]^S[33]^S[32]^S[17]^S[16]^S[1]^S[0]); break;
			case 50: S.set(position, S[48]^S[34]^S[32]^S[18]^S[16]^S[2]^S[0]); break;
			case 51: S.set(position, S[48]^S[34]^S[33]^S[18]^S[17]^S[3]^S[0]); break;
			case 52: S.set(position, S[48]^S[36]^S[32]^S[20]^S[16]^S[4]^S[0]); break;
			case 53: S.set(position, S[48]^S[36]^S[33]^S[20]^S[17]^S[5]^S[0]); break;
			case 54: S.set(position, S[48]^S[36]^S[34]^S[20]^S[18]^S[6]^S[0]); break;
			case 55: S.set(position, S[48]^S[36]^S[35]^S[20]^S[19]^S[7]^S[0]); break;
			case 56: S.set(position, S[48]^S[40]^S[32]^S[24]^S[16]^S[8]^S[0]); break;
			case 57: S.set(position, S[48]^S[40]^S[33]^S[24]^S[17]^S[9]^S[0]); break;
			case 58: S.set(position, S[48]^S[40]^S[34]^S[24]^S[18]^S[10]^S[0]); break;
			case 59: S.set(position, S[48]^S[40]^S[35]^S[24]^S[19]^S[11]^S[0]); break;
			case 60: S.set(position, S[48]^S[40]^S[36]^S[24]^S[20]^S[12]^S[0]); break;
			case 61: S.set(position, S[48]^S[40]^S[37]^S[24]^S[21]^S[13]^S[0]); break;
			case 62: S.set(position, S[48]^S[40]^S[38]^S[24]^S[22]^S[14]^S[0]); break;
			case 63: S.set(position, S[48]^S[40]^S[39]^S[24]^S[23]^S[15]^S[0]); break;
			case 67: S.set(position, S[66]^S[65]^S[64]^S[3]^S[2]^S[1]^S[0]); break;
			case 69: S.set(position, S[68]^S[65]^S[64]^S[5]^S[4]^S[1]^S[0]); break;
			case 70: S.set(position, S[68]^S[66]^S[64]^S[6]^S[4]^S[2]^S[0]); break;
			case 71: S.set(position, S[68]^S[66]^S[65]^S[6]^S[5]^S[3]^S[0]); break;
			case 73: S.set(position, S[72]^S[65]^S[64]^S[9]^S[8]^S[1]^S[0]); break;
			case 74: S.set(position, S[72]^S[66]^S[64]^S[10]^S[8]^S[2]^S[0]); break;
			case 75: S.set(position, S[72]^S[66]^S[65]^S[10]^S[9]^S[3]^S[0]); break;
			case 76: S.set(position, S[72]^S[68]^S[64]^S[12]^S[8]^S[4]^S[0]); break;
			case 77: S.set(position, S[72]^S[68]^S[65]^S[12]^S[9]^S[5]^S[0]); break;
			case 78: S.set(position, S[72]^S[68]^S[66]^S[12]^S[10]^S[6]^S[0]); break;
			case 79: S.set(position, S[72]^S[68]^S[67]^S[12]^S[11]^S[7]^S[0]); break;
			case 81: S.set(position, S[80]^S[65]^S[64]^S[17]^S[16]^S[1]^S[0]); break;
			case 82: S.set(position, S[80]^S[66]^S[64]^S[18]^S[16]^S[2]^S[0]); break;
			case 83: S.set(position, S[80]^S[66]^S[65]^S[18]^S[17]^S[3]^S[0]); break;
			case 84: S.set(position, S[80]^S[68]^S[64]^S[20]^S[16]^S[4]^S[0]); break;
			case 85: S.set(position, S[80]^S[68]^S[65]^S[20]^S[17]^S[5]^S[0]); break;
			case 86: S.set(position, S[80]^S[68]^S[66]^S[20]^S[18]^S[6]^S[0]); break;
			case 87: S.set(position, S[80]^S[68]^S[67]^S[20]^S[19]^S[7]^S[0]); break;
			case 88: S.set(position, S[80]^S[72]^S[64]^S[24]^S[16]^S[8]^S[0]); break;
			case 89: S.set(position, S[80]^S[72]^S[65]^S[24]^S[17]^S[9]^S[0]); break;
			case 90: S.set(position, S[80]^S[72]^S[66]^S[24]^S[18]^S[10]^S[0]); break;
			case 91: S.set(position, S[80]^S[72]^S[67]^S[24]^S[19]^S[11]^S[0]); break;
			case 92: S.set(position, S[80]^S[72]^S[68]^S[24]^S[20]^S[12]^S[0]); break;
			case 93: S.set(position, S[80]^S[72]^S[69]^S[24]^S[21]^S[13]^S[0]); break;
			case 94: S.set(position, S[80]^S[72]^S[70]^S[24]^S[22]^S[14]^S[0]); break;
			case 95: S.set(position, S[80]^S[72]^S[71]^S[24]^S[23]^S[15]^S[0]); break;
			case 97: S.set(position, S[96]^S[65]^S[64]^S[33]^S[32]^S[1]^S[0]); break;
			case 98: S.set(position, S[96]^S[66]^S[64]^S[34]^S[32]^S[2]^S[0]); break;
			case 99: S.set(position, S[96]^S[66]^S[65]^S[34]^S[33]^S[3]^S[0]); break;
			case 100: S.set(position, S[96]^S[68]^S[64]^S[36]^S[32]^S[4]^S[0]); break;
			case 101: S.set(position, S[96]^S[68]^S[65]^S[36]^S[33]^S[5]^S[0]); break;
			case 102: S.set(position, S[96]^S[68]^S[66]^S[36]^S[34]^S[6]^S[0]); break;
			case 103: S.set(position, S[96]^S[68]^S[67]^S[36]^S[35]^S[7]^S[0]); break;
			case 104: S.set(position, S[96]^S[72]^S[64]^S[40]^S[32]^S[8]^S[0]); break;
			case 105: S.set(position, S[96]^S[72]^S[65]^S[40]^S[33]^S[9]^S[0]); break;
			case 106: S.set(position, S[96]^S[72]^S[66]^S[40]^S[34]^S[10]^S[0]); break;
			case 107: S.set(position, S[96]^S[72]^S[67]^S[40]^S[35]^S[11]^S[0]); break;
			case 108: S.set(position, S[96]^S[72]^S[68]^S[40]^S[36]^S[12]^S[0]); break;
			case 109: S.set(position, S[96]^S[72]^S[69]^S[40]^S[37]^S[13]^S[0]); break;
			case 110: S.set(position, S[96]^S[72]^S[70]^S[40]^S[38]^S[14]^S[0]); break;
			case 111: S.set(position, S[96]^S[72]^S[71]^S[40]^S[39]^S[15]^S[0]); break;
			case 112: S.set(position, S[96]^S[80]^S[64]^S[48]^S[32]^S[16]^S[0]); break;
			case 113: S.set(position, S[96]^S[80]^S[65]^S[48]^S[33]^S[17]^S[0]); break;
			case 114: S.set(position, S[96]^S[80]^S[66]^S[48]^S[34]^S[18]^S[0]); break;
			case 115: S.set(position, S[96]^S[80]^S[67]^S[48]^S[35]^S[19]^S[0]); break;
			case 116: S.set(position, S[96]^S[80]^S[68]^S[48]^S[36]^S[20]^S[0]); break;
			case 117: S.set(position, S[96]^S[80]^S[69]^S[48]^S[37]^S[21]^S[0]); break;
			case 118: S.set(position, S[96]^S[80]^S[70]^S[48]^S[38]^S[22]^S[0]); break;
			case 119: S.set(position, S[96]^S[80]^S[71]^S[48]^S[39]^S[23]^S[0]); break;
			case 120: S.set(position, S[96]^S[80]^S[72]^S[48]^S[40]^S[24]^S[0]); break;
			case 121: S.set(position, S[96]^S[80]^S[73]^S[48]^S[41]^S[25]^S[0]); break;
			case 122: S.set(position, S[96]^S[80]^S[74]^S[48]^S[42]^S[26]^S[0]); break;
			case 123: S.set(position, S[96]^S[80]^S[75]^S[48]^S[43]^S[27]^S[0]); break;
			case 124: S.set(position, S[96]^S[80]^S[76]^S[48]^S[44]^S[28]^S[0]); break;
			case 125: S.set(position, S[96]^S[80]^S[77]^S[48]^S[45]^S[29]^S[0]); break;
			case 126: S.set(position, S[96]^S[80]^S[78]^S[48]^S[46]^S[30]^S[0]); break;
			case 127: S.set(position, S[96]^S[80]^S[79]^S[48]^S[47]^S[31]^S[0]); break;
			default: return false;
		}

	} else {
		return false;
	}

	return true;
}
#undef BITS_IN
#undef BITS_OUT

#if 0
#define BITS_IN 8
#define BITS_OUT 8
	 	// BITS=8 -- DEGREE 2
		// S[7] = S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
		// S[11] = S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0];
		// S[13] = S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0];
		// S[14] = S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0];
		// S[15] = S[12]^S[10]^S[9]^S[6]^S[5]^S[3]^S[0];
		// S[19] = S[18]^S[17]^S[16]^S[3]^S[2]^S[1]^S[0];
		// S[21] = S[20]^S[17]^S[16]^S[5]^S[4]^S[1]^S[0];
		// S[22] = S[20]^S[18]^S[16]^S[6]^S[4]^S[2]^S[0];
		// S[23] = S[20]^S[18]^S[17]^S[6]^S[5]^S[3]^S[0];
		// S[25] = S[24]^S[17]^S[16]^S[9]^S[8]^S[1]^S[0];
		// S[26] = S[24]^S[18]^S[16]^S[10]^S[8]^S[2]^S[0];
		// S[27] = S[24]^S[18]^S[17]^S[10]^S[9]^S[3]^S[0];
		// S[28] = S[24]^S[20]^S[16]^S[12]^S[8]^S[4]^S[0];
		// S[29] = S[24]^S[20]^S[17]^S[12]^S[9]^S[5]^S[0];
		// S[30] = S[24]^S[20]^S[18]^S[12]^S[10]^S[6]^S[0];
		// S[31] = S[24]^S[20]^S[19]^S[12]^S[11]^S[7]^S[0];
		// S[35] = S[34]^S[33]^S[32]^S[3]^S[2]^S[1]^S[0];
		// S[37] = S[36]^S[33]^S[32]^S[5]^S[4]^S[1]^S[0];
		// S[38] = S[36]^S[34]^S[32]^S[6]^S[4]^S[2]^S[0];
		// S[39] = S[36]^S[34]^S[33]^S[6]^S[5]^S[3]^S[0];
		// S[41] = S[40]^S[33]^S[32]^S[9]^S[8]^S[1]^S[0];
		// S[42] = S[40]^S[34]^S[32]^S[10]^S[8]^S[2]^S[0];
		// S[43] = S[40]^S[34]^S[33]^S[10]^S[9]^S[3]^S[0];
		// S[44] = S[40]^S[36]^S[32]^S[12]^S[8]^S[4]^S[0];
		// S[45] = S[40]^S[36]^S[33]^S[12]^S[9]^S[5]^S[0];
		// S[46] = S[40]^S[36]^S[34]^S[12]^S[10]^S[6]^S[0];
		// S[47] = S[40]^S[36]^S[35]^S[12]^S[11]^S[7]^S[0];
		// S[49] = S[48]^S[33]^S[32]^S[17]^S[16]^S[1]^S[0];
		// S[50] = S[48]^S[34]^S[32]^S[18]^S[16]^S[2]^S[0];
		// S[51] = S[48]^S[34]^S[33]^S[18]^S[17]^S[3]^S[0];
		// S[52] = S[48]^S[36]^S[32]^S[20]^S[16]^S[4]^S[0];
		// S[53] = S[48]^S[36]^S[33]^S[20]^S[17]^S[5]^S[0];
		// S[54] = S[48]^S[36]^S[34]^S[20]^S[18]^S[6]^S[0];
		// S[55] = S[48]^S[36]^S[35]^S[20]^S[19]^S[7]^S[0];
		// S[56] = S[48]^S[40]^S[32]^S[24]^S[16]^S[8]^S[0];
		// S[57] = S[48]^S[40]^S[33]^S[24]^S[17]^S[9]^S[0];
		// S[58] = S[48]^S[40]^S[34]^S[24]^S[18]^S[10]^S[0];
		// S[59] = S[48]^S[40]^S[35]^S[24]^S[19]^S[11]^S[0];
		// S[60] = S[48]^S[40]^S[36]^S[24]^S[20]^S[12]^S[0];
		// S[61] = S[48]^S[40]^S[37]^S[24]^S[21]^S[13]^S[0];
		// S[62] = S[48]^S[40]^S[38]^S[24]^S[22]^S[14]^S[0];
		// S[63] = S[48]^S[40]^S[39]^S[24]^S[23]^S[15]^S[0];
		// S[67] = S[66]^S[65]^S[64]^S[3]^S[2]^S[1]^S[0];
		// S[69] = S[68]^S[65]^S[64]^S[5]^S[4]^S[1]^S[0];
		// S[70] = S[68]^S[66]^S[64]^S[6]^S[4]^S[2]^S[0];
		// S[71] = S[68]^S[66]^S[65]^S[6]^S[5]^S[3]^S[0];
		// S[73] = S[72]^S[65]^S[64]^S[9]^S[8]^S[1]^S[0];
		// S[74] = S[72]^S[66]^S[64]^S[10]^S[8]^S[2]^S[0];
		// S[75] = S[72]^S[66]^S[65]^S[10]^S[9]^S[3]^S[0];
		// S[76] = S[72]^S[68]^S[64]^S[12]^S[8]^S[4]^S[0];
		// S[77] = S[72]^S[68]^S[65]^S[12]^S[9]^S[5]^S[0];
		// S[78] = S[72]^S[68]^S[66]^S[12]^S[10]^S[6]^S[0];
		// S[79] = S[72]^S[68]^S[67]^S[12]^S[11]^S[7]^S[0];
		// S[81] = S[80]^S[65]^S[64]^S[17]^S[16]^S[1]^S[0];
		// S[82] = S[80]^S[66]^S[64]^S[18]^S[16]^S[2]^S[0];
		// S[83] = S[80]^S[66]^S[65]^S[18]^S[17]^S[3]^S[0];
		// S[84] = S[80]^S[68]^S[64]^S[20]^S[16]^S[4]^S[0];
		// S[85] = S[80]^S[68]^S[65]^S[20]^S[17]^S[5]^S[0];
		// S[86] = S[80]^S[68]^S[66]^S[20]^S[18]^S[6]^S[0];
		// S[87] = S[80]^S[68]^S[67]^S[20]^S[19]^S[7]^S[0];
		// S[88] = S[80]^S[72]^S[64]^S[24]^S[16]^S[8]^S[0];
		// S[89] = S[80]^S[72]^S[65]^S[24]^S[17]^S[9]^S[0];
		// S[90] = S[80]^S[72]^S[66]^S[24]^S[18]^S[10]^S[0];
		// S[91] = S[80]^S[72]^S[67]^S[24]^S[19]^S[11]^S[0];
		// S[92] = S[80]^S[72]^S[68]^S[24]^S[20]^S[12]^S[0];
		// S[93] = S[80]^S[72]^S[69]^S[24]^S[21]^S[13]^S[0];
		// S[94] = S[80]^S[72]^S[70]^S[24]^S[22]^S[14]^S[0];
		// S[95] = S[80]^S[72]^S[71]^S[24]^S[23]^S[15]^S[0];
		// S[97] = S[96]^S[65]^S[64]^S[33]^S[32]^S[1]^S[0];
		// S[98] = S[96]^S[66]^S[64]^S[34]^S[32]^S[2]^S[0];
		// S[99] = S[96]^S[66]^S[65]^S[34]^S[33]^S[3]^S[0];
		// S[100] = S[96]^S[68]^S[64]^S[36]^S[32]^S[4]^S[0];
		// S[101] = S[96]^S[68]^S[65]^S[36]^S[33]^S[5]^S[0];
		// S[102] = S[96]^S[68]^S[66]^S[36]^S[34]^S[6]^S[0];
		// S[103] = S[96]^S[68]^S[67]^S[36]^S[35]^S[7]^S[0];
		// S[104] = S[96]^S[72]^S[64]^S[40]^S[32]^S[8]^S[0];
		// S[105] = S[96]^S[72]^S[65]^S[40]^S[33]^S[9]^S[0];
		// S[106] = S[96]^S[72]^S[66]^S[40]^S[34]^S[10]^S[0];
		// S[107] = S[96]^S[72]^S[67]^S[40]^S[35]^S[11]^S[0];
		// S[108] = S[96]^S[72]^S[68]^S[40]^S[36]^S[12]^S[0];
		// S[109] = S[96]^S[72]^S[69]^S[40]^S[37]^S[13]^S[0];
		// S[110] = S[96]^S[72]^S[70]^S[40]^S[38]^S[14]^S[0];
		// S[111] = S[96]^S[72]^S[71]^S[40]^S[39]^S[15]^S[0];
		// S[112] = S[96]^S[80]^S[64]^S[48]^S[32]^S[16]^S[0];
		// S[113] = S[96]^S[80]^S[65]^S[48]^S[33]^S[17]^S[0];
		// S[114] = S[96]^S[80]^S[66]^S[48]^S[34]^S[18]^S[0];
		// S[115] = S[96]^S[80]^S[67]^S[48]^S[35]^S[19]^S[0];
		// S[116] = S[96]^S[80]^S[68]^S[48]^S[36]^S[20]^S[0];
		// S[117] = S[96]^S[80]^S[69]^S[48]^S[37]^S[21]^S[0];
		// S[118] = S[96]^S[80]^S[70]^S[48]^S[38]^S[22]^S[0];
		// S[119] = S[96]^S[80]^S[71]^S[48]^S[39]^S[23]^S[0];
		// S[120] = S[96]^S[80]^S[72]^S[48]^S[40]^S[24]^S[0];
		// S[121] = S[96]^S[80]^S[73]^S[48]^S[41]^S[25]^S[0];
		// S[122] = S[96]^S[80]^S[74]^S[48]^S[42]^S[26]^S[0];
		// S[123] = S[96]^S[80]^S[75]^S[48]^S[43]^S[27]^S[0];
		// S[124] = S[96]^S[80]^S[76]^S[48]^S[44]^S[28]^S[0];
		// S[125] = S[96]^S[80]^S[77]^S[48]^S[45]^S[29]^S[0];
		// S[126] = S[96]^S[80]^S[78]^S[48]^S[46]^S[30]^S[0];
		// S[127] = S[96]^S[80]^S[79]^S[48]^S[47]^S[31]^S[0];
		// S[131] = S[130]^S[129]^S[128]^S[3]^S[2]^S[1]^S[0];
		// S[133] = S[132]^S[129]^S[128]^S[5]^S[4]^S[1]^S[0];
		// S[134] = S[132]^S[130]^S[128]^S[6]^S[4]^S[2]^S[0];
		// S[135] = S[132]^S[130]^S[129]^S[6]^S[5]^S[3]^S[0];
		// S[137] = S[136]^S[129]^S[128]^S[9]^S[8]^S[1]^S[0];
		// S[138] = S[136]^S[130]^S[128]^S[10]^S[8]^S[2]^S[0];
		// S[139] = S[136]^S[130]^S[129]^S[10]^S[9]^S[3]^S[0];
		// S[140] = S[136]^S[132]^S[128]^S[12]^S[8]^S[4]^S[0];
		// S[141] = S[136]^S[132]^S[129]^S[12]^S[9]^S[5]^S[0];
		// S[142] = S[136]^S[132]^S[130]^S[12]^S[10]^S[6]^S[0];
		// S[143] = S[136]^S[132]^S[131]^S[12]^S[11]^S[7]^S[0];
		// S[145] = S[144]^S[129]^S[128]^S[17]^S[16]^S[1]^S[0];
		// S[146] = S[144]^S[130]^S[128]^S[18]^S[16]^S[2]^S[0];
		// S[147] = S[144]^S[130]^S[129]^S[18]^S[17]^S[3]^S[0];
		// S[148] = S[144]^S[132]^S[128]^S[20]^S[16]^S[4]^S[0];
		// S[149] = S[144]^S[132]^S[129]^S[20]^S[17]^S[5]^S[0];
		// S[150] = S[144]^S[132]^S[130]^S[20]^S[18]^S[6]^S[0];
		// S[151] = S[144]^S[132]^S[131]^S[20]^S[19]^S[7]^S[0];
		// S[152] = S[144]^S[136]^S[128]^S[24]^S[16]^S[8]^S[0];
		// S[153] = S[144]^S[136]^S[129]^S[24]^S[17]^S[9]^S[0];
		// S[154] = S[144]^S[136]^S[130]^S[24]^S[18]^S[10]^S[0];
		// S[155] = S[144]^S[136]^S[131]^S[24]^S[19]^S[11]^S[0];
		// S[156] = S[144]^S[136]^S[132]^S[24]^S[20]^S[12]^S[0];
		// S[157] = S[144]^S[136]^S[133]^S[24]^S[21]^S[13]^S[0];
		// S[158] = S[144]^S[136]^S[134]^S[24]^S[22]^S[14]^S[0];
		// S[159] = S[144]^S[136]^S[135]^S[24]^S[23]^S[15]^S[0];
		// S[161] = S[160]^S[129]^S[128]^S[33]^S[32]^S[1]^S[0];
		// S[162] = S[160]^S[130]^S[128]^S[34]^S[32]^S[2]^S[0];
		// S[163] = S[160]^S[130]^S[129]^S[34]^S[33]^S[3]^S[0];
		// S[164] = S[160]^S[132]^S[128]^S[36]^S[32]^S[4]^S[0];
		// S[165] = S[160]^S[132]^S[129]^S[36]^S[33]^S[5]^S[0];
		// S[166] = S[160]^S[132]^S[130]^S[36]^S[34]^S[6]^S[0];
		// S[167] = S[160]^S[132]^S[131]^S[36]^S[35]^S[7]^S[0];
		// S[168] = S[160]^S[136]^S[128]^S[40]^S[32]^S[8]^S[0];
		// S[169] = S[160]^S[136]^S[129]^S[40]^S[33]^S[9]^S[0];
		// S[170] = S[160]^S[136]^S[130]^S[40]^S[34]^S[10]^S[0];
		// S[171] = S[160]^S[136]^S[131]^S[40]^S[35]^S[11]^S[0];
		// S[172] = S[160]^S[136]^S[132]^S[40]^S[36]^S[12]^S[0];
		// S[173] = S[160]^S[136]^S[133]^S[40]^S[37]^S[13]^S[0];
		// S[174] = S[160]^S[136]^S[134]^S[40]^S[38]^S[14]^S[0];
		// S[175] = S[160]^S[136]^S[135]^S[40]^S[39]^S[15]^S[0];
		// S[176] = S[160]^S[144]^S[128]^S[48]^S[32]^S[16]^S[0];
		// S[177] = S[160]^S[144]^S[129]^S[48]^S[33]^S[17]^S[0];
		// S[178] = S[160]^S[144]^S[130]^S[48]^S[34]^S[18]^S[0];
		// S[179] = S[160]^S[144]^S[131]^S[48]^S[35]^S[19]^S[0];
		// S[180] = S[160]^S[144]^S[132]^S[48]^S[36]^S[20]^S[0];
		// S[181] = S[160]^S[144]^S[133]^S[48]^S[37]^S[21]^S[0];
		// S[182] = S[160]^S[144]^S[134]^S[48]^S[38]^S[22]^S[0];
		// S[183] = S[160]^S[144]^S[135]^S[48]^S[39]^S[23]^S[0];
		// S[184] = S[160]^S[144]^S[136]^S[48]^S[40]^S[24]^S[0];
		// S[185] = S[160]^S[144]^S[137]^S[48]^S[41]^S[25]^S[0];
		// S[186] = S[160]^S[144]^S[138]^S[48]^S[42]^S[26]^S[0];
		// S[187] = S[160]^S[144]^S[139]^S[48]^S[43]^S[27]^S[0];
		// S[188] = S[160]^S[144]^S[140]^S[48]^S[44]^S[28]^S[0];
		// S[189] = S[160]^S[144]^S[141]^S[48]^S[45]^S[29]^S[0];
		// S[190] = S[160]^S[144]^S[142]^S[48]^S[46]^S[30]^S[0];
		// S[191] = S[160]^S[144]^S[143]^S[48]^S[47]^S[31]^S[0];
		// S[193] = S[192]^S[129]^S[128]^S[65]^S[64]^S[1]^S[0];
		// S[194] = S[192]^S[130]^S[128]^S[66]^S[64]^S[2]^S[0];
		// S[195] = S[192]^S[130]^S[129]^S[66]^S[65]^S[3]^S[0];
		// S[196] = S[192]^S[132]^S[128]^S[68]^S[64]^S[4]^S[0];
		// S[197] = S[192]^S[132]^S[129]^S[68]^S[65]^S[5]^S[0];
		// S[198] = S[192]^S[132]^S[130]^S[68]^S[66]^S[6]^S[0];
		// S[199] = S[192]^S[132]^S[131]^S[68]^S[67]^S[7]^S[0];
		// S[200] = S[192]^S[136]^S[128]^S[72]^S[64]^S[8]^S[0];
		// S[201] = S[192]^S[136]^S[129]^S[72]^S[65]^S[9]^S[0];
		// S[202] = S[192]^S[136]^S[130]^S[72]^S[66]^S[10]^S[0];
		// S[203] = S[192]^S[136]^S[131]^S[72]^S[67]^S[11]^S[0];
		// S[204] = S[192]^S[136]^S[132]^S[72]^S[68]^S[12]^S[0];
		// S[205] = S[192]^S[136]^S[133]^S[72]^S[69]^S[13]^S[0];
		// S[206] = S[192]^S[136]^S[134]^S[72]^S[70]^S[14]^S[0];
		// S[207] = S[192]^S[136]^S[135]^S[72]^S[71]^S[15]^S[0];
		// S[208] = S[192]^S[144]^S[128]^S[80]^S[64]^S[16]^S[0];
		// S[209] = S[192]^S[144]^S[129]^S[80]^S[65]^S[17]^S[0];
		// S[210] = S[192]^S[144]^S[130]^S[80]^S[66]^S[18]^S[0];
		// S[211] = S[192]^S[144]^S[131]^S[80]^S[67]^S[19]^S[0];
		// S[212] = S[192]^S[144]^S[132]^S[80]^S[68]^S[20]^S[0];
		// S[213] = S[192]^S[144]^S[133]^S[80]^S[69]^S[21]^S[0];
		// S[214] = S[192]^S[144]^S[134]^S[80]^S[70]^S[22]^S[0];
		// S[215] = S[192]^S[144]^S[135]^S[80]^S[71]^S[23]^S[0];
		// S[216] = S[192]^S[144]^S[136]^S[80]^S[72]^S[24]^S[0];
		// S[217] = S[192]^S[144]^S[137]^S[80]^S[73]^S[25]^S[0];
		// S[218] = S[192]^S[144]^S[138]^S[80]^S[74]^S[26]^S[0];
		// S[219] = S[192]^S[144]^S[139]^S[80]^S[75]^S[27]^S[0];
		// S[220] = S[192]^S[144]^S[140]^S[80]^S[76]^S[28]^S[0];
		// S[221] = S[192]^S[144]^S[141]^S[80]^S[77]^S[29]^S[0];
		// S[222] = S[192]^S[144]^S[142]^S[80]^S[78]^S[30]^S[0];
		// S[223] = S[192]^S[144]^S[143]^S[80]^S[79]^S[31]^S[0];
		// S[224] = S[192]^S[160]^S[128]^S[96]^S[64]^S[32]^S[0];
		// S[225] = S[192]^S[160]^S[129]^S[96]^S[65]^S[33]^S[0];
		// S[226] = S[192]^S[160]^S[130]^S[96]^S[66]^S[34]^S[0];
		// S[227] = S[192]^S[160]^S[131]^S[96]^S[67]^S[35]^S[0];
		// S[228] = S[192]^S[160]^S[132]^S[96]^S[68]^S[36]^S[0];
		// S[229] = S[192]^S[160]^S[133]^S[96]^S[69]^S[37]^S[0];
		// S[230] = S[192]^S[160]^S[134]^S[96]^S[70]^S[38]^S[0];
		// S[231] = S[192]^S[160]^S[135]^S[96]^S[71]^S[39]^S[0];
		// S[232] = S[192]^S[160]^S[136]^S[96]^S[72]^S[40]^S[0];
		// S[233] = S[192]^S[160]^S[137]^S[96]^S[73]^S[41]^S[0];
		// S[234] = S[192]^S[160]^S[138]^S[96]^S[74]^S[42]^S[0];
		// S[235] = S[192]^S[160]^S[139]^S[96]^S[75]^S[43]^S[0];
		// S[236] = S[192]^S[160]^S[140]^S[96]^S[76]^S[44]^S[0];
		// S[237] = S[192]^S[160]^S[141]^S[96]^S[77]^S[45]^S[0];
		// S[238] = S[192]^S[160]^S[142]^S[96]^S[78]^S[46]^S[0];
		// S[239] = S[192]^S[160]^S[143]^S[96]^S[79]^S[47]^S[0];
		// S[240] = S[192]^S[160]^S[144]^S[96]^S[80]^S[48]^S[0];
		// S[241] = S[192]^S[160]^S[145]^S[96]^S[81]^S[49]^S[0];
		// S[242] = S[192]^S[160]^S[146]^S[96]^S[82]^S[50]^S[0];
		// S[243] = S[192]^S[160]^S[147]^S[96]^S[83]^S[51]^S[0];
		// S[244] = S[192]^S[160]^S[148]^S[96]^S[84]^S[52]^S[0];
		// S[245] = S[192]^S[160]^S[149]^S[96]^S[85]^S[53]^S[0];
		// S[246] = S[192]^S[160]^S[150]^S[96]^S[86]^S[54]^S[0];
		// S[247] = S[192]^S[160]^S[151]^S[96]^S[87]^S[55]^S[0];
		// S[248] = S[192]^S[160]^S[152]^S[96]^S[88]^S[56]^S[0];
		// S[249] = S[192]^S[160]^S[153]^S[96]^S[89]^S[57]^S[0];
		// S[250] = S[192]^S[160]^S[154]^S[96]^S[90]^S[58]^S[0];
		// S[251] = S[192]^S[160]^S[155]^S[96]^S[91]^S[59]^S[0];
		// S[252] = S[192]^S[160]^S[156]^S[96]^S[92]^S[60]^S[0];
		// S[253] = S[192]^S[160]^S[157]^S[96]^S[93]^S[61]^S[0];
		// S[254] = S[192]^S[160]^S[158]^S[96]^S[94]^S[62]^S[0];
		// S[255] = S[192]^S[160]^S[159]^S[96]^S[95]^S[63]^S[0];
#undef BITS_IN
#undef BITS_OUT
#endif
} /* namespace apn */


/*
 * Return the degree of the N polynomials
 */
template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::compute_degree() {
	const function_in_out_t& S = *this;

	vector<int> degrees;

	for (size_t i = 0; i < constants_out_t::WORD_BITS; i++) {
	 
	 	bitset_in_t polynomials;

		for (word_in_t j = 0; ; j++) {
			if ((S[j]>>i) & 1){
				for (word_in_t k = j; ; k++) {
					if ((j&k) == j) {
						polynomials.flip(k);
					}
					if (k == constants_in_t::WORD_MAX)
						break;
				}/*k*/
			}
			if (j == constants_in_t::WORD_MAX)
				break;
		}/*j*/

		int degree = 0;
		for(word_in_t j = 0; ; j++) {
			int weight_j = __builtin_popcount(j);
			if (polynomials[j] && (degree<weight_j)) {
				degree = weight_j;
			}
			if (j == constants_in_t::WORD_MAX)
				break;
		}/*j*/

		degrees.push_back(degree);
	
	}/*i*/

	degree = *max_element(degrees.begin(), degrees.end());

	return degree;
}

/*
 * Return the degrees of monomials of the Sbox
 */
/*
template <size_t BITS_IN, size_t BITS_OUT>
std::vector<int> function_in_out_t::get_degrees() const {
	return vector<int>();
	//return degrees;
}
*/

/*
 * Return the total degree of the Sbox
 */
template <size_t BITS_IN, size_t BITS_OUT>
int function_in_out_t::get_degree() const {
	return degree;
}

#endif

