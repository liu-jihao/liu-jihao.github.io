// Text input and output of GMP rationals for viprchk.cpp: parses and prints integers and
// fractions with mpq_set_str and mpq_get_str, rejects malformed tokens and zero denominators.
// It does not change the checking logic.
#include <iostream>
#include <string>
#include <vector>
#include <gmp.h>
std::istream& operator>>(std::istream& in, mpq_ptr q) {
  std::string s; if(!(in>>s)) return in;
  if(mpq_set_str(q,s.c_str(),10)!=0 || mpz_sgn(mpq_denref(q))==0) {
    in.setstate(std::ios::failbit); return in;
  }
  mpq_canonicalize(q); return in;
}
std::ostream& operator<<(std::ostream& out, mpq_srcptr q) {
  std::vector<char> s(mpz_sizeinbase(mpq_numref(q),10)+mpz_sizeinbase(mpq_denref(q),10)+4);
  mpq_get_str(s.data(),10,q); return out<<s.data();
}
