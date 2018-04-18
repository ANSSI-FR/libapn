#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <algorithm>
#include <bitset>
#include <vector>
#include <map>
#include <iostream>
#include <string>

using namespace std;

#define N        256
#define DEGREE   3

// uint8_t prod(uint8_t A[6], uint8_t x) {
// 	uint8_t val = 0;
// 	for(int i=0; i<6; ++i) {
// 		val <<= 1;
// 		val |= __builtin_popcount(A[i]&x)&1;
// 	}/*i*/
// 	return val;
// }

// void test() {

// 	uint8_t S0[64] = {
// 		0x00,0x36,0x30,0x0d,0x0f,0x12,0x35,0x23,0x19,0x3f,0x2d,0x34,0x03,0x14,0x29,0x21,
// 		0x3b,0x24,0x02,0x22,0x0a,0x08,0x39,0x25,0x3c,0x13,0x2a,0x0e,0x32,0x1a,0x3a,0x18,
// 		0x27,0x1b,0x15,0x11,0x10,0x1d,0x01,0x3e,0x2f,0x28,0x33,0x38,0x07,0x2b,0x2c,0x26,
// 		0x1f,0x0b,0x04,0x1c,0x3d,0x2e,0x05,0x31,0x09,0x06,0x17,0x20,0x1e,0x0c,0x37,0x16,
// 	};

// 	uint8_t S1[64];

// 	//S0(x) = B(S1(A(x)+9)+4
// 	//B^-1(S0(x)+4) = S1(A(x)+9)

// 	uint8_t A[6] = {
// 		0b110101,
// 		0b111100,
// 		0b100000,
// 		0b000101,
// 		0b000100,
// 		0b000110,
// 	};

// 	// uint8_t B[6][6] = {
// 	// 	{0, 1, 1, 1, 0, 1},
// 	// 	{0, 0, 0, 0, 0, 1},
// 	// 	{0, 0, 1, 1, 1, 0},
// 	// 	{0, 0, 0, 1, 1, 1},
// 	// 	{0, 0, 1, 0, 1, 0},
// 	// 	{1, 0, 1, 1, 0, 1},
// 	// };

// 	uint8_t invB[6] = {
// 		0b000111,
// 		0b100110,
// 		0b011100,
// 		0b001010,
// 		0b011110,
// 		0b010000,
// 	};


// 	uint8_t a = 9;
// 	uint8_t b = 4;

// 	//B^-1(S0(x)+4) = S1(A(x)+9)
// 	for(int x=0; x<64; ++x) {
// 		uint8_t Ax9  = prod(A, x)^9;
// 		uint8_t S0x4 = prod(invB, S0[x]^4);

// 		S1[Ax9] = S0x4;
// 	}/*x*/

// 	printf("s = mq.SBox([");
// 	for(int i=0; i<64; ++i) printf("%d,", S1[i]);
// 	printf("])\n");
// }







int main() {

	/* 6-bit involutive APN permutation affine-equivalent to Dillon's one */
	// uint8_t S1={0,19,14,61,47,36,49,26,33,58,11,10,56,39,2,15,29,20,40,1,17,41,22,30,62,31,7,60,28,16,23,25,55,8,37,35,5,34,51,13,18,21,52,48,46,45,44,4,43,6,50,38,42,54,53,32,12,57,9,63,27,3,24,59};

	// test(); return 0;

	vector< bitset<N> > eq;

#if DEGREE==2
	for(uint64_t a=1; a<N; ++a) {
		for(uint64_t b=a+1; b<N; ++b) {

			for(uint64_t x=0; x<N; ++x) {
				for(uint64_t y=x+1; y<N; ++y) {

					if(x>y) continue;
					if(!(x<(x^a))) continue;
					if(!(x<(x^b))) continue;
					if(!(x<(x^a^b))) continue;

					if(!(y<(y^a))) continue;
					if(!(y<(y^b))) continue;
					if(!(y<(y^a^b))) continue;

					bitset<N> s;
	
					s.set(x);
					s.set(x^a);
					s.set(x^b);
					s.set(x^a^b);

					s.set(y);
					s.set(y^a);
					s.set(y^b);
					s.set(y^a^b);

					eq.push_back(s);
				}
			}
		}
	}

#elif DEGREE==3

	for(uint64_t a=0; a<N; ++a) {
		for(uint64_t b=a+1; b<N; ++b) {
			for(uint64_t c=b+1; c<N; ++c) {

				for(uint64_t x=0; x<N; ++x) {
					for(uint64_t y=x+1; y<N; ++y) {

						if(x>y) continue;

						if(!(x<(x^a))) continue;
						if(!(x<(x^b))) continue;
						if(!(x<(x^c))) continue;
						if(!(x<(x^a^b))) continue;
						if(!(x<(x^a^c))) continue;
						if(!(x<(x^b^c))) continue;
						if(!(x<(x^a^b^c))) continue;

						if(!(y<(y^a))) continue;
						if(!(y<(y^b))) continue;
						if(!(y<(y^c))) continue;
						if(!(y<(y^a^b))) continue;
						if(!(y<(y^a^c))) continue;
						if(!(y<(y^b^c))) continue;
						if(!(y<(y^a^b^c))) continue;

						bitset<N> s;

						s.set(x);
						s.set(x^a);
						s.set(x^b);
						s.set(x^c);
						s.set(x^a^b);
						s.set(x^a^c);
						s.set(x^b^c);
						s.set(x^a^b^c);

						s.set(y);
						s.set(y^a);
						s.set(y^b);
						s.set(y^c);
						s.set(y^a^b);
						s.set(y^a^c);
						s.set(y^b^c);
						s.set(y^a^b^c);

						eq.push_back(s);
					}
					
				}
			}
		}
	}

#elif DEGREE==4

	for(uint64_t a=1; a<N; ++a) {
		for(uint64_t b=a+1; b<N; ++b) {
			for(uint64_t c=b+1; c<N; ++c) {
				for(uint64_t d=c+1; d<N; ++d) {

					for(uint64_t x=0; x<N; ++x) {
						for(uint64_t y=x+1; y<N; ++y) {

							if(x>=y) continue;

							if(!(x<(x^a))) continue;
							if(!(x<(x^b))) continue;
							if(!(x<(x^c))) continue;
							if(!(x<(x^d))) continue;
							if(!(x<(x^a^b))) continue;
							if(!(x<(x^a^c))) continue;
							if(!(x<(x^a^d))) continue;
							if(!(x<(x^b^c))) continue;
							if(!(x<(x^b^d))) continue;
							if(!(x<(x^c^d))) continue;
							if(!(x<(x^a^b^c))) continue;
							if(!(x<(x^a^b^d))) continue;
							if(!(x<(x^a^c^d))) continue;
							if(!(x<(x^b^c^d))) continue;
							if(!(x<(x^a^b^c^d))) continue;

							if(!(y<(y^a))) continue;
							if(!(y<(y^b))) continue;
							if(!(y<(y^c))) continue;
							if(!(y<(y^d))) continue;
							if(!(y<(y^a^b))) continue;
							if(!(y<(y^a^c))) continue;
							if(!(y<(y^a^d))) continue;
							if(!(y<(y^b^c))) continue;
							if(!(y<(y^b^d))) continue;
							if(!(y<(y^c^d))) continue;
							if(!(y<(y^a^b^c))) continue;
							if(!(y<(y^a^b^d))) continue;
							if(!(y<(y^a^c^d))) continue;
							if(!(y<(y^b^c^d))) continue;
							if(!(y<(y^a^b^c^d))) continue;

							bitset<N> s;

							s.set(x);
							s.set(x^a);
							s.set(x^b);
							s.set(x^c);
							s.set(x^d);
							s.set(x^a^b);
							s.set(x^a^c);
							s.set(x^a^d);
							s.set(x^b^c);
							s.set(x^b^d);
							s.set(x^c^d);
							s.set(x^a^b^c);
							s.set(x^a^b^d);
							s.set(x^a^c^d);
							s.set(x^b^c^d);
							s.set(x^a^b^c^d);

							s.set(y);
							s.set(y^a);
							s.set(y^b);
							s.set(y^c);
							s.set(y^d);
							s.set(y^a^b);
							s.set(y^a^c);
							s.set(y^a^d);
							s.set(y^b^c);
							s.set(y^b^d);
							s.set(y^c^d);
							s.set(y^a^b^c);
							s.set(y^a^b^d);
							s.set(y^a^c^d);
							s.set(y^b^c^d);
							s.set(y^a^b^c^d);

							eq.push_back(s);
						}
					}
				}
			}
		}
	}

#endif

	sort(eq.begin(), eq.end(), 
		[](const bitset<N>& lhs, const bitset<N>& rhs) {
			return lhs.to_string() < rhs.to_string();
		});

	eq.erase(
		unique(eq.begin(), eq.end()),
		eq.end()
	);

	map<int, bitset<N> > selected;

	for(auto x: eq) {
		int leading = N;
		for(int b=N-1; b>=0; --b) {
			if(x.test(b)) {
				if(leading==N) {
					leading=b;
					if(selected.find(leading)==selected.end()) {
						selected[leading] = x;
						break;
					}
				}
			}
		}
	}

	cout << "// N=" << N << " -- DEGREE=" << DEGREE << endl;

	for(auto it=selected.begin();
	 it!=selected.end(); ++it) {


		cout << "S[" << it->first << "] = ";

	 	bool first = true;

	 	std::string s ("");
		for(int b=N-1; b>=0; --b) {
			if(it->second.test(b)) {

				if(first)  first=false;
				else       s += "S[" + to_string(b) + "]^";

			}
		}/*b*/

		s.pop_back();
		cout << s+";" << endl;
	}

	// // 5-bit degree-2 permutations
	// S[ 7] = S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
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

	// // 5-bit degree-3 permutations
	// S[15] = S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[23] = S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[27] = S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0];
	// S[29] = S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0];
	// S[30] = S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0];
	// S[31] = S[28]^S[26]^S[25]^S[22]^S[21]^S[19]^S[16]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1];



	// // 6-bit degree-2 permutations
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



	// // 6-bit degree-3 permutations
	// S[15] = S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[23] = S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[27] = S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0];
	// S[29] = S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0];
	// S[30] = S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0];
	// S[31] = S[28]^S[26]^S[25]^S[22]^S[21]^S[19]^S[16]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1];
	// S[39] = S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[43] = S[42]^S[41]^S[40]^S[35]^S[34]^S[33]^S[32]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0];
	// S[45] = S[44]^S[41]^S[40]^S[37]^S[36]^S[33]^S[32]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0];
	// S[46] = S[44]^S[42]^S[40]^S[38]^S[36]^S[34]^S[32]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0];
	// S[47] = S[44]^S[42]^S[41]^S[38]^S[37]^S[35]^S[32]^S[14]^S[13]^S[11]^S[8]^S[7]^S[4]^S[2]^S[1];
	// S[51] = S[50]^S[49]^S[48]^S[35]^S[34]^S[33]^S[32]^S[19]^S[18]^S[17]^S[16]^S[3]^S[2]^S[1]^S[0];
	// S[53] = S[52]^S[49]^S[48]^S[37]^S[36]^S[33]^S[32]^S[21]^S[20]^S[17]^S[16]^S[5]^S[4]^S[1]^S[0];
	// S[54] = S[52]^S[50]^S[48]^S[38]^S[36]^S[34]^S[32]^S[22]^S[20]^S[18]^S[16]^S[6]^S[4]^S[2]^S[0];
	// S[55] = S[52]^S[50]^S[49]^S[38]^S[37]^S[35]^S[32]^S[22]^S[21]^S[19]^S[16]^S[7]^S[4]^S[2]^S[1];
	// S[57] = S[56]^S[49]^S[48]^S[41]^S[40]^S[33]^S[32]^S[25]^S[24]^S[17]^S[16]^S[9]^S[8]^S[1]^S[0];
	// S[58] = S[56]^S[50]^S[48]^S[42]^S[40]^S[34]^S[32]^S[26]^S[24]^S[18]^S[16]^S[10]^S[8]^S[2]^S[0];
	// S[59] = S[56]^S[50]^S[49]^S[42]^S[41]^S[35]^S[32]^S[26]^S[25]^S[19]^S[16]^S[11]^S[8]^S[2]^S[1];
	// S[60] = S[56]^S[52]^S[48]^S[44]^S[40]^S[36]^S[32]^S[28]^S[24]^S[20]^S[16]^S[12]^S[8]^S[4]^S[0];
	// S[61] = S[56]^S[52]^S[49]^S[44]^S[41]^S[37]^S[32]^S[28]^S[25]^S[21]^S[16]^S[13]^S[8]^S[4]^S[1];
	// S[62] = S[56]^S[52]^S[50]^S[44]^S[42]^S[38]^S[32]^S[28]^S[26]^S[22]^S[16]^S[14]^S[8]^S[4]^S[2];
	// S[63] = S[56]^S[52]^S[51]^S[44]^S[43]^S[39]^S[32]^S[28]^S[27]^S[23]^S[16]^S[15]^S[8]^S[4]^S[3];



	// // 6-bit degree-4 permutations
	// S[31] = S[30]^S[29]^S[28]^S[27]^S[26]^S[25]^S[24]^S[23]^S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[15]^S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[47] = S[46]^S[45]^S[44]^S[43]^S[42]^S[41]^S[40]^S[39]^S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[15]^S[14]^S[13]^S[12]^S[11]^S[10]^S[9]^S[8]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[55] = S[54]^S[53]^S[52]^S[51]^S[50]^S[49]^S[48]^S[39]^S[38]^S[37]^S[36]^S[35]^S[34]^S[33]^S[32]^S[23]^S[22]^S[21]^S[20]^S[19]^S[18]^S[17]^S[16]^S[7]^S[6]^S[5]^S[4]^S[3]^S[2]^S[1]^S[0];
	// S[59] = S[58]^S[57]^S[56]^S[51]^S[50]^S[49]^S[48]^S[43]^S[42]^S[41]^S[40]^S[35]^S[34]^S[33]^S[32]^S[27]^S[26]^S[25]^S[24]^S[19]^S[18]^S[17]^S[16]^S[11]^S[10]^S[9]^S[8]^S[3]^S[2]^S[1]^S[0];
	// S[61] = S[60]^S[57]^S[56]^S[53]^S[52]^S[49]^S[48]^S[45]^S[44]^S[41]^S[40]^S[37]^S[36]^S[33]^S[32]^S[29]^S[28]^S[25]^S[24]^S[21]^S[20]^S[17]^S[16]^S[13]^S[12]^S[9]^S[8]^S[5]^S[4]^S[1]^S[0];
	// S[62] = S[60]^S[58]^S[56]^S[54]^S[52]^S[50]^S[48]^S[46]^S[44]^S[42]^S[40]^S[38]^S[36]^S[34]^S[32]^S[30]^S[28]^S[26]^S[24]^S[22]^S[20]^S[18]^S[16]^S[14]^S[12]^S[10]^S[8]^S[6]^S[4]^S[2]^S[0];
	// S[63] = S[60]^S[58]^S[57]^S[54]^S[53]^S[51]^S[48]^S[46]^S[45]^S[43]^S[40]^S[39]^S[36]^S[34]^S[33]^S[30]^S[29]^S[27]^S[24]^S[23]^S[20]^S[18]^S[17]^S[15]^S[12]^S[10]^S[9]^S[6]^S[5]^S[3]^S[0];

	return 0;
}
