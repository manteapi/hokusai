#ifndef ALJABR_MATRIX_HPP
#define ALJABR_MATRIX_HPP

#include <array>
#include <algorithm>
#include <iostream>

namespace aljabr
{
    ///Class representing a matrix with N rows and M columns of a type T
    template <class T, size_t N, size_t M>
        class Matrix
        {
            public :

                std::array<T, N*M> d; ///< Data container

                //! Default constructor.
                /*!
                 * \brief By default all elements are set to zero.
                 */
                Matrix(){ d.fill( T(0.0) ); }

                //! Constructor with initial value.
                /*!
                 * \brief Constructor using a user defined default value to fill the container.
                 * \param value default value
                 */
                Matrix(const T& value){ d.fill(value); }

                //! Clone constructor
                Matrix(const Matrix& m){ d = m.array(); }

                //! Default destructor.
                ~Matrix(){}

                //! Return the size of the container.
                /*!
                 * \brief Return the size of the container.
                 * \return The size of the container.
                 */
                size_t size() const { return N*M;}

                //! Return the number of rows.
                /*!
                 * \brief Return the number of rows.
                 * \return The number of rows.
                 */
                size_t rows() const { return N;}

                //! Return the number of columns.
                /*!
                 * \brief Return the number of columns.
                 * \return The number of columns.
                 */
                size_t cols() const { return M;}

                //! Return a pointer to the first element of the data container.
                /*!
                 * \brief Return a pointer to the first element of the data container.
                 * \return A pointer to the first element of the data container.
                 */
                T* data(){ return d.data(); }

                //! Return a const pointer to the first element of the data container.
                /*!
                 * \brief Return a const pointer to the first element of the data container.
                 * \return A const pointer to the first element of the data container.
                 */
                const T* data() const { return d.data(); }

                //! Return a reference to the container
                std::array<T, N*M> & array() { return d; }

                //! Return a const reference to the container
                const std::array<T, N*M> & array() const { return d; }

                //! Return the L2 norm of the vector.
                /*!
                 * \brief Return the L2 norm of the vector.
                 * \return L2 norm of the vector.
                 */
                T norm2()
                {
                    T l2 = 0;
                    for(size_t i=0; i<size(); ++i)
                        l2 += d[i]*d[i];
                    return sqrt(l2);
                }

                //! Normalize value of the container using L2 norm.
                /*!
                 * \brief Normalize value of the container using L2 norm.
                 */
                void normalize()
                {
                    T l2 = norm2();
                    for(size_t i=0; i<size(); ++i)
                        d[i]/=l2;
                }

                //! Accessor to the ith element of the container.
                /*!
                 * \brief Accessor to the ith element of the container.
                 * \return A reference to the ith element.
                 * \param i index.
                 */
                T& operator[](size_t i){ return d[i]; }
                
                //! Const accessor to the ith element of the container.
                /*!
                 * \brief Const accessor to the ith element of the container.
                 * \return A const reference to the ith element.
                 * \param i index.
                 */
                const T& operator[](size_t i) const { return d[i]; }
                
                //! Accessor to the (ith, jth) element of the matrix
                /*!
                 * \brief Accessor to the (ith, jth) element of the matrix
                 * \return A reference to the (ith, jth) element.
                 * \param i line index
                 * \param j column index
                 */
                T& operator()(size_t i, size_t j){ return d[i*M+j]; }

                //! Const accessor to the (ith, jth) element of the matrix
                /*!
                 * \brief Const accessor to the (ith, jth) element of the matrix
                 * \return A const reference to the (ith, jth) element.
                 * \param i line index
                 * \param j column index
                 */
                const T& operator()(size_t i, size_t j) const { return d[i*M+j]; }

                //! Build an identity matrix
                void identity()
                {
                    for(size_t i=0; i<N; ++i)
                    {
                        for(size_t j=0; j<M; ++j)
                        {
                            if(i==j)
                                d[i*M+j] = 1.0;
                            else
                                d[i*M+j] = 0.0;
                        }
                    }
                }

                //! Operator equal.
                /*!
                 * \brief Copy a matrix using the operator equal.
                 * \param m matrix to copy.
                 * \return A reference to itself.
                 */
                Matrix& operator=(const Matrix& m){ d = m.array(); return *this;}

                //! Output stream function.
                /*!
                 * \brief Utility function to display matrix's content.
                 */
                friend std::ostream& operator<<( std::ostream& os, const Matrix& m)
                {
                    for(size_t i=0; i<N; ++i)
                    {
                        for(size_t j=0; j<M; ++j)
                        {
                            os << m(i,j) << " ";
                        }
                        os << std::endl;
                    }
                    return os;
                }
        };
}//namespace aljabr

#endif //ALJABR_MATRIX_HPP
