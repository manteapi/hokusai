#ifndef ALJABR_UTILS_INL
#define ALJABR_UTILS_INL

namespace aljabr
{

template<class T>
int signum(T val)
{
    return (T(0)<val)-(val<T(0));
}

}//namespace aljabr

#endif //ALJABR_UTILS_INL
