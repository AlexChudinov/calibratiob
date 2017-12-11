#include <Eigen/Dense>
#include <algorithm>
#include "regression.h"

using namespace Eigen;
using namespace std;

Regression::Regression(const XYWData &xywData, size_t degree)
    :
      m_bOk(true)
{
    if(xywData.size() <= degree)
    {
        m_bOk = false;
    }
    else
    {
        m_vCoefs.assign(degree + 1, 0.0);
        calculateCoefs(xywData);
    }
}

Regression::Regression()
    :
      m_bOk(true)
{

}

bool Regression::isValid() const
{
    return m_bOk;
}

Regression::operator bool() const
{
    return isValid();
}

double Regression::yValue(double fXValue) const
{
    double ret = 0.0;
    double fXPowN = 1.0;
    for(double c : m_vCoefs)
    {
        ret += c * fXPowN;
        fXPowN *= fXValue;
    }
    return ret;
}

const Regression::Vector& Regression::coefs() const
{
    return m_vCoefs;
}

void Regression::coefs(const Vector &c)
{
    m_vCoefs = c;
}

void Regression::valid(bool bIsValid)
{
    m_bOk = bIsValid;
}

void Regression::calculateCoefs(const XYWData &xywData)
{
    MatrixXd M(xywData.size(), m_vCoefs.size()), W;
    VectorXd Y(xywData.size());

    size_t i = 0;
    W.setZero(xywData.size(), xywData.size());
    for(const auto& xyw : xywData)
    {
        Y(i)  = xyw.second.first;
        W(i,i)= xyw.second.second;
        double dXPowN = 1.0, dXVal = xyw.first;
        for(size_t j = 0; j < m_vCoefs.size(); ++j)
        {
            M(i, j) = dXPowN;
            dXPowN *= dXVal;
        }
        ++i;
    }

    MatrixXd TM = M.transpose() * W;
    MatrixXd vCoefs = (TM * M).ldlt().solve(TM * Y);

    copy(vCoefs.data(), vCoefs.data() + vCoefs.rows(), m_vCoefs.begin());
}

LinearRegression::LinearRegression(const XYWData &xywData, const VectorFun &vFuns)
    :
      Regression()
{
    if(xywData.size() < vFuns.size())
    {
        valid(false);
    }
    else
    {
        m_vFuns = vFuns;
        calculateCoefs(xywData);
    }
}

LinearRegression::LinearRegression
(
    const Vector &x,
    const Vector &y,
    const Vector &w,
    const VectorFun &vF
)
    :
      m_vFuns(vF)
{
    if(vF.size() > x.size())
    {
        valid(false);
    }
    else
    {
        assert(x.size() == y.size() && x.size() == w.size());
        MatrixXd M(x.size(), vF.size()), W;
        VectorXd Y(x.size());
        W.setZero(x.size(), x.size());
        for(size_t i = 0; i < x.size(); ++i)
        {
            Y(i) = y[i];
            W(i, i) = w[i];
            for(size_t j = 0; j < vF.size(); ++j)
            {
                M(i, j) = vF[j](x[i]);
            }
        }
        MatrixXd TM = M.transpose() * W;
        MatrixXd vCoefs = (TM*M).ldlt().solve(TM*Y);
        coefs(Vector(vCoefs.data(), vCoefs.data() + vCoefs.rows()));
    }
}

double LinearRegression::yValue(double fXValue) const
{
    double ret = 0;
    for(size_t i = 0; i < coefs().size(); ++i)
    {
        ret += coefs()[i] * m_vFuns[i](fXValue);
    }
    return ret;
}

void LinearRegression::calculateCoefs(const XYWData &xywData)
{
    MatrixXd M(xywData.size(), m_vFuns.size()), W;
    VectorXd Y(xywData.size());

    size_t i = 0;
    W.setZero(xywData.size(), xywData.size());
    for(const auto& xyw : xywData)
    {
        Y(i)  = xyw.second.first;
        W(i,i)= xyw.second.second;
        for(size_t j = 0; j < m_vFuns.size(); ++j)
        {
            M(i, j) = m_vFuns[j](xyw.first);
        }
        ++i;
    }

    MatrixXd TM = M.transpose() * W;
    MatrixXd vCoefs = (TM * M).ldlt().solve(TM * Y);

    coefs(Vector(vCoefs.data(), vCoefs.data() + vCoefs.rows()));
}
