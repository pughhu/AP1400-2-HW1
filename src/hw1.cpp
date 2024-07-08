#include "hw1.h"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

Matrix algebra::zeros(size_t n, size_t m){
    // return Matrix{};
    std::vector<std::vector<double>> m1(n,std::vector<double>(m,0));
    return m1;
}

Matrix algebra::ones(size_t n, size_t m){
    std::vector<std::vector<double>> m1(n,std::vector<double>(m,1));
    return m1;
}

Matrix algebra::random(size_t n, size_t m, double min, double max){
    if(min > max)   throw std::logic_error("min cannot be greater than max!");
    std::default_random_engine dre;
    std::uniform_real_distribution<double> di(min, max);
    
    // std::vector<std::vector<double>> m1(n,std::vector<double>(m,di(dre)));
    
    std::vector<std::vector<double>> m1(n,std::vector<double>(m));
    // for(std::size_t i=0; i<n; ++i){
    //     for(std::size_t j=0; j<m; ++j){
    //         m1[i][j] = di(dre);
    //     }
    // }
    for(auto& e1: m1){
        for(auto& e2: e1){
            e2 = di(dre);
        }
    }
    return m1;
}


void algebra::show(const Matrix& matrix){
    int precision = 3;
    for(auto vec: matrix){
        for(auto num: vec){
            std::cout<<std::fixed<<std::setprecision(precision)<<num<<"  ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

Matrix algebra::multiply(const Matrix& matrix, double c){
    Matrix m1(matrix);
    // const std::size_t m = matrix.size();
    // const std::size_t n = matrix.at(0).size();
    
    // for(std::size_t i=0; i<m; ++i){
    //     for(std::size_t j=0; j<n; ++j){
    //         m1[i][j] *= c;
    //     }
    // }

    for(auto& e1: m1){
        for(auto& e2: e1){
            e2 *= c;
        }
    }

    return m1;
}

Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2){
    if(matrix1.empty() || matrix2.empty())  return Matrix{};

    if(matrix1[0].size() != matrix2.size()) throw std::logic_error("dimension is not identical!");

    const size_t m = matrix1.size();
    const size_t n = matrix1[0].size();
    const size_t p = matrix2[0].size();

    std::vector<std::vector<double>> m1(m,std::vector<double>(p,0));

    for(size_t i=0; i<m; ++i){
        for(size_t j=0; j<p; ++j){
            for(size_t k=0; k<n; ++k){
                m1[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }

    return m1;
}

Matrix algebra::sum(const Matrix& matrix, double c){
    if(matrix.empty())  return Matrix{};

    Matrix m1(matrix);
    // const auto m = m1.size();
    // const auto n = m1[0].size();

    // for(auto i=0; i<m; ++i){
    //     for(auto j=0; j<n; ++j){
    //         m1[i][j] += c;
    //     }
    // }

    for(auto& e1:m1){
        for(auto& e2:e1){
            e2 += c;
        }
    }

    return m1;
}

Matrix algebra::sum(const Matrix& matrix1, const Matrix& matrix2){
    if(matrix1.empty() && matrix2.empty())  
        return Matrix{};
    else if(matrix1.empty() || matrix2.empty()) 
        throw std::logic_error{"dimension fault"};
    const auto m1 = matrix1.size();
    const auto n1 = matrix1[0].size();
    const auto m2 = matrix1.size();
    const auto n2 = matrix2[0].size();

    if(m1 != m2 || n1 != n2)    throw std::logic_error("dimension is not identical!");

    Matrix mt(m1,std::vector<double>(n1,0));
    for(auto i=0; i<m1; ++i){
        for(auto j=0; j<n1; ++j){
            mt[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    return mt;
}

Matrix algebra::transpose(const Matrix& matrix){
    if(matrix.empty())  return Matrix{};
    
    const auto m = matrix.size();
    const auto n = matrix[0].size();

    Matrix mt(n,std::vector<double>(m,0));
    for(auto i=0; i<n; ++i){
        for(auto j=0; j<m; ++j){
            mt[i][j] = matrix[j][i];
        }
    }

    return mt;
}

Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m){
    if(matrix.empty())  throw std::logic_error("no dimension");
    auto n1 = matrix.size();
    auto m1 = matrix[0].size();
    if(n>=n1 || m>=m1)  throw std::logic_error("dimension fault");
    if(n1==1 || m1==1)  return Matrix{};

    auto lamba = [n,m](size_t i, size_t j){
        size_t row,col;
        if(i<n)
            row = i;
        else
            row = i+1;
        if(j<m)
            col = j;
        else
            col = j+1;
        return std::pair<size_t,size_t>(row,col);
    };

    Matrix mt(n1-1,std::vector<double>(m1-1,0));
    for(auto i=0; i<n1-1; ++i){
        for(auto j=0; j<m1-1; ++j){
            // auto [row, col] = lamba(i,j);
            auto resutl = lamba(i,j);
            auto row = resutl.first;
            auto col = resutl.second;
            mt[i][j] = matrix[row][col];
        }
    }

    return mt;
}

double algebra::determinant(const Matrix& matrix){
    if(matrix.empty())  return 1;
    if(matrix.size() != matrix[0].size())   throw std::logic_error("dimension not identical");

    auto lamba = [](size_t i, size_t j){
        if((i+j)%2 == 0)
            return 1;
        else
            return -1;
    };

    auto n = matrix.size();
    auto m = matrix[0].size();
    double result{};
    for(auto i=0; i<n; ++i){
        result += lamba(i,0)*matrix[i][0]*algebra::determinant(algebra::minor(matrix,i,0));
    }
    
    return result;
}

Matrix algebra::inverse(const Matrix& matrix){
    if(matrix.empty())   return Matrix{};
    const size_t n = matrix.size();
    const size_t m = matrix[0].size();
    if(n != m)  throw std::logic_error("matrix n is not equal to m");

    // debug
    // algebra::show(matrix);
    // construct a identical matrix
    Matrix I(n,std::vector<double>(n,0));
    for(auto i=0; i<n; ++i){
        I[i][i] = 1;
    }

    Matrix merge(n,std::vector<double>(2*n,0));
    for(auto i=0; i<n; ++i){
        for(auto j=0; j<2*n; ++j){
            if(j<n)
                merge[i][j] = matrix[i][j];
            else
                merge[i][j] = I[i][j-n];
        }
    }

    // debug
    // algebra::show(merge);
    // guass eliminatation
    auto select_main = [&](int step){
        size_t i = step;
        while(i<merge.size()){
            if(merge[i][step] != 0){
                break;
            }else{
                ++i;
            }
        }
        if(i == merge.size())  throw std::logic_error("no inverse.");
        if(i != step){
            // mt[i].swap(mt.[step]);
            swap(merge[i], merge[step]);
        }

        auto old = merge[step][step];        
        for(auto i=step; i<merge[0].size(); ++i){
            merge[step][i] /= old;
            // algebra::show(merge);
        }

        //
        for(auto i=step+1; i<merge.size(); ++i){
            auto k = merge[i][step];
            for(auto j=step; j<merge[0].size(); ++j){
                merge[i][j] -= k * merge[step][j];
            }
        }
    };

    // iterate until all done
    for(int i=0; i<n; ++i){
        select_main(i);
        // algebra::show(merge);
    }

    // iterate reverse
    for(int col=n-1; col>0; --col){
        for(int row=col-1; row>=0; --row){
            auto k = merge[row][col];
            for(int i=0; i<merge[0].size(); ++i){
                merge[row][i] -= k * merge[col][i];
            }
        }
        // algebra::show(merge);
    }

    // finished
    Matrix result(n, std::vector<double>(n,0));
    for(auto i=0; i<n; ++i){
        for(auto j=0; j<n; ++j){
            result[i][j] = merge[i][j+n];
        }
    }

    return result;
}

Matrix algebra::concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis){
    if(matrix1.empty()) return matrix2;
    if(matrix2.empty()) return matrix1;

    auto n1 = matrix1.size();
    auto m1 = matrix1[0].size();
    auto n2 = matrix2.size();
    auto m2 = matrix2[0].size();
    if(axis == 0){
        if(m1 != m2)
            throw std::logic_error("dimension not idential");
        Matrix mt(n1+n2, std::vector<double>(m2,0));
        for(auto i=0; i<n1+n2; ++i){
            mt[i] = i<n1 ? matrix1[i] : matrix2[i-n1];
        }

        return mt;
    }else{
        if(n1 != n2)
            throw std::logic_error("dimension not identical");
        // Matrix mt(n1,std::vector<double>(m1+m2,0));
        Matrix mt(n1,std::vector<double>());
        for(auto i=0; i<n1; ++i){
            mt[i].reserve(m1+m2);
        }
        for(auto i=0; i<n1; ++i){
            // mt[i].assign(matrix1[i].begin(), matrix1[i].end());
            mt[i].insert(mt[i].begin(), matrix1[i].begin(), matrix1[i].end());
            mt[i].insert(mt[i].end(), matrix2[i].begin(), matrix2[i].end());
        }

        return mt;
    }
}

Matrix algebra::ero_swap(const Matrix& matrix, size_t r1, size_t r2){
    auto size = matrix.size();
    auto max = [](size_t r1, size_t r2){return r1>r2 ? r1: r2;};
    if(max(r1,r2) >= size)  throw std::logic_error("size error");

    Matrix mt(matrix);
    // swap(mt[r1],mt[r2]);
    mt[r1].swap(mt[r2]);
    return mt;
}

Matrix algebra::ero_multiply(const Matrix& matrix, size_t r, double c){
    if(r >= matrix.size())  throw std::logic_error("size error");
    
    Matrix mt(matrix);
    for(auto i=0; i<mt[r].size(); ++i){
        mt[r][i] *= c;
    }

    return mt;
}

Matrix algebra::ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2){
    auto size = matrix.size();
    auto max = [](size_t r1, size_t r2){return r1>r2 ? r1: r2;};
    if(max(r1,r2) >= size)  throw std::logic_error("size error");

    Matrix mt(matrix);
    for(auto i=0; i<mt[r2].size(); ++i){
        mt[r2][i] += c* mt[r1][i];
    }

    return mt;
}

Matrix algebra::upper_triangular(const Matrix& matrix){
    if(matrix.empty())  return Matrix{};
    if(matrix.size() != matrix[0].size())   throw std::logic_error("size error");
    Matrix mt(matrix);
    // swap
    auto n = mt.size();
    auto m = mt[0].size();
    for(auto i=0; i<n; ++i){
        auto main = i;
        while(main<n && mt[main][i] == 0){
            ++main;
        }
        mt = algebra::ero_swap(mt, i, main);
        for(auto j=i+1; j<n; ++j){
            auto k = (mt[j][i] / mt[i][i]) * (-1);
            mt = algebra::ero_sum(mt,i,k,j);
        }
    }

    return mt;
}