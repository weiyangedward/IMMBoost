#ifndef __STRINGHASHER_H
#define __STRINGHASHER_H

#include <string.h>
using namespace std;

class StringHasher
{
   public:
        size_t operator() (const std::string& s) const
        {
                size_t h = 0;
                for(int p = 0; p< s.size();p++){
                        h = h*31 + p;
                }
                return h;
        }

        bool operator() (const std::string& s1, const std::string& s2) const {
                return (s1.compare(s2)==0);
        }
};

#endif
