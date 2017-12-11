#ifndef REGRESSION_H
#define REGRESSION_H

#include <map>
#include <vector>
#include <functional>

/**
 * @brief The Regression class implements polynomial regression
 */
class Regression
{
public:
    typedef std::vector<double> Vector;
    typedef std::pair<double, double> WeightedYVal;
    typedef std::pair<double, WeightedYVal> XYWPoint;
    typedef std::map<double, WeightedYVal> XYWData;
    typedef std::size_t size_t;

    Regression(const XYWData& xywData, size_t degree);

    const Vector& coefs() const;

    void coefs(const Vector& c);

    bool isValid() const;

    /**
     * @brief operator bool
     * @return true if Regression is valid
     */
    operator bool() const;

    /**
     * @brief yValue
     * @param fXValue
     * @return yValue calculated using regression coefficient
     */
    double yValue(double fXValue) const;

protected:
    Regression();

    void valid(bool bIsValid = true);

private:
    Vector m_vCoefs;
    bool m_bOk;

    void calculateCoefs(const XYWData& xywData);
};

/**
 * @brief The LinearRegression class implements linear regression of any type
 */
class LinearRegression : public Regression
{
public:
    typedef std::vector< std::function<double(double)> > VectorFun;

    LinearRegression(const XYWData& xywData, const VectorFun& vFuns);

    /**
     * @brief LinearRegression calculates regression from raw data values
     * @param x
     * @param y
     * @param w statistic weights
     * @param vF vector of regression basic functions
     */
    LinearRegression
    (
        const Vector& x,
        const Vector& y,
        const Vector& w,
        const VectorFun& vF
    );

    double yValue(double fXValue) const;

private:
    VectorFun m_vFuns;

    void calculateCoefs(const XYWData& xywData);
};

#endif // REGRESSION_H
