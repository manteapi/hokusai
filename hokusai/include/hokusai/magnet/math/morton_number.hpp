#pragma once

#include <stdint.h>
#include <limits>
#include "dilated_int.hpp"

namespace magnet {
    namespace math {
        /*! \brief A class for calculating Morton numbers.
         *
         * This class contains d DilatedIntegers, and uses them to perform
         * Morton Ordered addressing.
         *
         * \tparam d The dilation (or dimensionality) of the MortonNumber.
         */
        template<size_t d>
            struct MortonNumber
            {
                //! \brief Default constructor.
                inline MortonNumber() {}

                //! \brief Construct a MortonNumber from an actual morton integer
                inline MortonNumber(const size_t& MortonNum)
                {
                    for (size_t i(0); i < d; ++i)
                        _data[i].setDilatedValue(MortonNum >> i);
                }


                //! \brief Helper constructor for 3D MortonNumber's.
                inline MortonNumber(const size_t& x, const size_t& y, const size_t& z)
                {
                    _data[0] = x;
                    _data[1] = y;
                    _data[2] = z;
                }

                //! \brief Helper constructor for 3D MortonNumber's
                inline MortonNumber(const DilatedInteger<d>& x,
                        const DilatedInteger<d>& y,
                        const DilatedInteger<d>& z)
                {
                    _data[0] = x;
                    _data[1] = y;
                    _data[2] = z;
                }

                //! \brief Returns the morton number stored in this class.
                inline size_t getMortonNum()
                {
                    size_t retval = _data[0].getDilatedValue();
                    for (size_t i(1); i < d; ++i)
                        retval += _data[i].getDilatedValue() << i;

                    return retval;
                }

                /*! \brief Accessor to the component DilatedInteger classes
                 * forming the MortonNumber.
                 */
                inline const DilatedInteger<d>& operator[](size_t i) const
                {
                    return _data[i];
                }

                /*! \brief Accessor to the component DilatedInteger classes
                 * forming the MortonNumber.
                 */
                inline DilatedInteger<d>& operator[](size_t i)
                {
                    return _data[i];
                }

                /*! \brief Simple math operation for adding two morton numbers.
                 */
                inline MortonNumber operator+(const MortonNumber& o) const
                {
                    MortonNumber retval;
                    for (size_t i(0); i < d; ++i)
                        retval._data[i] = o._data[i] + _data[i];

                    return retval;
                }

                protected:
                DilatedInteger<d> _data[d];
            };
    }
}
