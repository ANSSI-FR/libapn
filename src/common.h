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

#ifndef COMMON_H
#define COMMON_H

#include <cstring>
#include <cstdint>
#include <climits>
#include <cinttypes>
#include <cmath>

#include <type_traits>

#include <array>
#include <bitset>

#if defined(BITS) || defined(BITS_IN) || defined(BITS_OUT)
#error "BITS* shoudl not be defined before including common.h"
#endif

// algorithms tweaking stuff
#define SKIP_CLEANUP 1

// common.h internal stuff
#define small_uint_types uint8_t, uint16_t, uint32_t, uint64_t
#define min_uint_type uint8_t
#define max_uint_type uint64_t
#define TWOTOTHE(A) (1<<(A))

// convinience external stuff
#define ffsull FFSULL<BITS>
#define ffs FFS<BITS>
#define ffs_in FFS<BITS_IN>
#define ffs_out FFS<BITS_OUT>
#define print_word PRINT_WORD<BITS>
#define print_word_in PRINT_WORD<BITS_IN>
#define print_word_out PRINT_WORD<BITS_OUT>
#define print_hword PRINT_HWORD<BITS>
#define print_hword_in PRINT_HWORD<BITS_IN>
#define print_hword_out PRINT_HWORD<BITS_OUT>
#define word_t word<BITS>
#define word_in_t word<BITS_IN>
#define word_out_t word<BITS_OUT>
#define xor_t xword<BITS>
#define xor_in_t xword<BITS_IN>
#define xor_out_t xword<BITS_OUT>
#define constants_t constants<BITS>
#define constants_in_t constants<BITS_IN>
#define constants_out_t constants<BITS_OUT>
#define bitset_t wbitset<BITS>
#define bitset_in_t wbitset<BITS_IN>
#define bitset_out_t wbitset<BITS_OUT>
#define array_t warray<BITS, BITS>
#define array_in_out_t warray<BITS_IN, BITS_OUT>
#define bbt_t wbbt<BITS_IN, BITS_OUT>
#define bbt_in_out_t wbbt<BITS_IN, BITS_OUT>
#define ddt_t wddt<BITS_IN, BITS_OUT>
#define ddt_in_out_t wddt<BITS_IN, BITS_OUT>
#define function_t function<BITS, BITS>
#define function_in_out_t function<BITS_IN, BITS_OUT>
#define linear_function_t linear_function<BITS, BITS>
#define linear_function_in_out_t linear_function<BITS_IN, BITS_OUT>
#define permutation_t permutation<BITS>
#define permutation_in_t permutation<BITS_IN>
#define permutation_out_t permutation<BITS_OUT>
#define linear_permutation_t linear_permutation<BITS>
#define linear_permutation_in_t linear_permutation<BITS_IN>
#define linear_permutation_out_t linear_permutation<BITS_OUT>
#define involution_t involution<BITS>
#define involution_in_t involution<BITS_IN>
#define involution_out_t involution<BITS_OUT>

using namespace::std;

namespace apn {

	/*
	// Easy solution
	// Type for words
	template <int N>
	using word = typename std::conditional<N <= 8, uint8_t, uint16_t>::type;
	// Type for xors
	template <int N>
	using xword = typename std::conditional<N == 3, uint32_t, typename std::conditional<N == 4, uint64_t, __uint128_t>::type>::type;
	*/

	// More general solution
	//C++14
	template <typename T>
	struct wrap_type {
		using type = T;
	};
	// Base case
	template <size_t BITS, bool fail, typename... Ts>
	struct smallest_uint_type : wrap_type<max_uint_type> {static_assert(!fail || BITS <= sizeof(max_uint_type)*CHAR_BIT, "The number of bits is too large.");};
	// Recursive case
	template <size_t BITS, bool fail, typename T, typename... Ts>
	struct smallest_uint_type<BITS, fail, T, Ts...> : std::conditional_t<BITS <= sizeof(T)*CHAR_BIT, wrap_type<T>, smallest_uint_type<BITS, fail, Ts...>> {};
	// Type for words
	template <size_t BITS>
	using word = typename smallest_uint_type<BITS, true, small_uint_types>::type;
	// Type for xors
	template <size_t BITS>
	using xword = typename smallest_uint_type<sizeof(word_t)*CHAR_BIT*TWOTOTHE(BITS)/2, false, small_uint_types>::type;

	template <size_t BITS>
	struct constants {
		static constexpr word_t WORD_MAX = (word_t) (TWOTOTHE(BITS)-1);//std::numeric_limits<word_t>::max(); // numeric_limits<word_t>::digits(2);
		static constexpr size_t WORD_BITS = sizeof(word_t)*CHAR_BIT; // numeric_limits<word_t>::digits(2);
		static constexpr size_t XOR_BITS = sizeof(xor_t)*CHAR_BIT;
		static constexpr size_t xCHARS = (WORD_BITS + 3) / 4; // std::ceil(std::log(N)/std::log(16));
		static constexpr size_t XCHARS = xCHARS;
		static constexpr size_t uCHARS = std::numeric_limits<word_t>::digits10 + 1; // std::ceil(std::log(N)/std::log(10) + 1);
		static constexpr size_t dCHARS = uCHARS;
		static constexpr size_t WORD_CHARS = uCHARS;
		static constexpr size_t HWORD_CHARS = XCHARS + 2;
		//static constexpr char WORD_FORMAT[] = sprintf();
		//static constexpr char HWORD_FORMAT[] = sprintf();
	};

	template<size_t BITS>
	using wbitset = typename std::bitset<TWOTOTHE(BITS)>;


 	template <size_t BITS>
	int FFSULL(const bitset_t& x) {
		return __builtin_ffsll((x).to_ullong())-1;
	}

 	template <size_t BITS>
	int FFS(const bitset_t& x) {
 		const bitset_t mask = ~(0ULL);
		bitset_t xull = x & mask;
		int res = 0;
		size_t shift = 0;

		for (; shift < constants_t::WORD_BITS; shift += sizeof(unsigned long long int)) {
			xull = (x >> shift) & mask;
			if (xull.any()) {
				res = shift + ffsull(xull);
				break;
			}
		}

		return res;
	}

	template <size_t BITS>
	void PRINT_WORD(char* buffer, word_t x) {
		sprintf(buffer, "%*" PRIuMAX, (int) constants_t::uCHARS, (uintmax_t) x);
	}

	template <size_t BITS>
	void PRINT_HWORD(char* buffer, word_t x) {
		sprintf(buffer, "0x%0*" PRIXMAX, (int) constants_t::XCHARS, (uintmax_t) x);
	}

	template <size_t BITS_IN, size_t BITS_OUT>
	using warray = typename std::array<word_out_t, TWOTOTHE(BITS_IN)>;

	template <size_t BITS_IN, size_t BITS_OUT>
	using wbbt = typename std::array<bitset_out_t, TWOTOTHE(BITS_IN)>;

	template <size_t BITS_IN, size_t BITS_OUT>
	using wddt = typename std::array<array_in_out_t, TWOTOTHE(BITS_IN)>;
}

#endif
