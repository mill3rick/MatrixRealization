//
// Created by mille on 29.02.2020.
//

#ifndef UNTITLED20_MAIN_HPP
#define UNTITLED20_MAIN_HPP
using namespace std;
template <typename t>
class Matrix;
template <typename t>
ostream& operator << (ostream& strm,Matrix<t> c);
template <typename t>
istream& operator >> (istream& strm,Matrix<t> &c);
template <typename t>
class Matrix
{
public:
    vector<t> row;
    vector<vector<t>> col;
    Matrix(size_t rows, size_t cols);
    Matrix() = default;
    Matrix<t> operator +(Matrix<t> b);
    Matrix<t> operator *(Matrix<t> b);
    Matrix<t> operator *(t b);
    void inverse();
    void transpose();
    void RowEchelonForm();
    void ReduceRowEchelonForm();
    void OMN();
    double det();
    friend ostream& operator <<<t> (ostream& strm,Matrix c);
    friend istream& operator >><t> (istream& strm,Matrix &c);
};

template <typename t>
Matrix<t>::Matrix(size_t rows, size_t cols)
{
    row.resize(rows);
    for(size_t i=0;i<cols;i++)
        col.push_back(row);
}

template <typename t>
Matrix<t> Matrix<t>::operator +(Matrix<t> b)
{
    Matrix c(row.size(),col.size());
    if(col.size()!=b.col.size() || col.at(0).size()!=b.col.at(0).size())
    {
        throw "ERR: Matrichi raznih razmerov";

    }

    for(size_t i = 0; i < col.size();i++)
    {
        for(size_t z = 0;z < row.size();z++)
        {
            c.col.at(i).at(z) = b.col.at(i).at(z)+col.at(i).at(z);
        }
    }
    c.OMN();
    return c;
}

template <typename t>
Matrix<t> Matrix<t>::operator *(t b)
{
    Matrix c(row.size(),col.size());
    for(size_t i = 0; i < col.size();i++)
    {
        for(size_t z = 0;z < row.size();z++)
        {
            c.col.at(i).at(z) = col.at(i).at(z) * b;
        }
    }
    c.OMN();
    return c;
}

template <typename t>
void Matrix<t>::transpose()
{
    Matrix<t> c(col.size(),row.size());
    for(size_t i = 0; i < col.size();i++)
    {
        for(size_t k = 0; k < row.size();k++)
            c.col.at(k).at(i) = col.at(i).at(k);
    }
    row.resize(c.col.at(0).size());
    for(size_t i=0;i<c.col.size();i++)
        col.push_back(row);
    col = c.col;
    OMN();
}
template <typename t>
ostream& operator << (ostream& strm,Matrix<t> c)
{
    c.OMN();
    for(size_t i = 0; i < c.col.size();i++)
    {
        for(size_t z = 0;z<c.row.size();z++)
        {
            strm << c.col.at(i).at(z)<<" ";
        }
        strm << endl;
    }
    return strm;
}
template <typename t>
istream& operator >> (istream& strm,Matrix<t> &c)
{
    for(size_t i = 0; i < c.col.size();i++)
    {
        for(size_t z = 0;z<c.row.size();z++)
        {
            strm >> c.col.at(i).at(z);
        }
    }
    return strm;
}

template <typename t>
Matrix<t> Matrix<t>::operator *(Matrix<t> b)
{
    Matrix c(row.size(),col.size());
    if(col.size()!=b.col.at(0).size())
    {
        throw "ERR: Matrichi raznih razmerov";

    }

    for (size_t i = 0; i < col.size(); i++)
    {
        for (size_t j = 0; j < row.size(); j++)
        {
            c.col[i][j] = 0;
            for (size_t q = 0; q < row.size(); q++)
            {
                c.col[i][j] += b.col[i][q] * col[q][j];
            }
        }
    }
    c.OMN();
    return c;
}
template <typename t>
double Matrix<t>::det()
{
    auto a = col;
    for(size_t i = 0; i < row.size() - 1; i++)
        for(size_t k = i + 1; k < row.size(); k++)
        {
            double coeff = -a[k][i] / a[i][i];
            for(size_t z = i; z < row.size(); z++)
                a[k][z] += a[i][z] * coeff;
        }
    double Det = 1;
    for(size_t i = 0; i < row.size(); i++)
        Det *= a[i][i];
    if( to_string(Det)=="nan" || to_string(Det)=="-nan" || Det==-0)
        Det=0;

    return Det;
}
template <typename t>
void Matrix<t>::RowEchelonForm()
{

    t q;
    for (size_t i = 0; i < col.size(); i++)
    {
        for (size_t z = i; z < col.at(0).size(); z++)
        {	if (col.at(z).at(i) != 0)
            {
                swap(col.at(i), col.at(z));
                break;
            }
        }

        q = col.at(i).at(i);

        for (size_t z = 0; z < col.at(0).size(); z++)
        {
            if (i < col.at(0).size() && q != 0) col.at(i).at(z) /= q;
        }

        for (size_t z = 0; z < col.size(); z++)
        {
            q = col.at(z).at(i);
            for (size_t k = 0; k < col.at(0).size(); k++)
            {
                if (z > i)
                    col.at(z).at(k) -= col.at(i).at(k) * q;
            }

        }
    }
    OMN();
    return ;
}
template <typename t>
void Matrix<t>::ReduceRowEchelonForm()
{
    t q;
    RowEchelonForm();
    for(size_t i = 0; i < col.size()-1;i++)
    {
        q = col.at(i).at(i+1);
        for(size_t k = 0; k < col.at(0).size();k++)
            col.at(i).at(k) = col.at(i).at(k) - col.at(i+1).at(k)*q;

    }
    OMN();
    return;
}
template <typename t>
void Matrix<t>::inverse()
{

    if(det()==0)
    {
        throw "Err: Det=0";
    }

    Matrix<t> l(col.at(0).size(),col.size());
    for(size_t i=0;i<l.col.size();i++)
    {
        l.col.at(i).at(i) = 1;
    }

    t q;//ql;
    for (size_t i = 0; i < col.size(); i++)
    {
        for (size_t z = i; z < col.at(0).size(); z++)
        {	if (col.at(z).at(i) != 0)
            {
                swap(col.at(i), col.at(z));
                swap(l.col.at(i), l.col.at(z));
                break;
            }
        }

        q = col.at(i).at(i);
        for (size_t z = 0; z < col.at(0).size(); z++)
        {
            if (i < col.at(0).size() && q != 0)
            {
                col.at(i).at(z) /= q;
                l.col.at(i).at(z) /= q;
            }
        }

        for (size_t z = 0; z < col.size(); z++)
        {
            q = col.at(z).at(i);
            //ql = l.col.at(z).at(i);
            for (size_t k = 0; k < col.at(0).size(); k++)
            {
                if (z > i)
                {
                    col.at(z).at(k) -= col.at(i).at(k) * q;
                    l.col.at(z).at(k) -= l.col.at(i).at(k) * q;
                }
            }

        }
    }
    for(size_t i = 0; i < col.size()-1;i++)
    {
        q = col.at(i).at(i+1);
        for(size_t k = 0; k < col.at(0).size();k++)
        {
            col.at(i).at(k) = col.at(i).at(k) - col.at(i+1).at(k)*q;
            l.col.at(i).at(k) = l.col.at(i).at(k) - l.col.at(i+1).at(k)*q;
        }
    }

    col = l.col;
    OMN();
    return;
}


template <typename t>
void Matrix<t>::OMN()
{
    for(size_t i = 0; i < col.size();i++)
    {
        for(size_t z = 0;z < row.size();z++)
        {
            if(col.at(i).at(z)==-0)
                col.at(i).at(z) = 0;
        }
    }

    return;
}


#endif //UNTITLED20_MAIN_HPP
