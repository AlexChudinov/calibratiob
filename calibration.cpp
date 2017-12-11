#include <cmath>
#include <cassert>
#include <algorithm>

#include "calibration.h"


Calibration::CalTypesList Calibration::s_supportedTypes
{
    {"Regular two point", RegularTwoPoint},
    {"2nd order polynomial", SecondOrderPoly},
    {"3rd order polynomial", ThirdOrderPoly},
    {"4th order 3 point", ForthDegreeThreePoint},
    {"bth order 3 point", ExpThreePoint}
};

Calibration* Calibration::create(CalibrationType type, const MassTimeIntensity &table)
{
    switch(type)
    {
    case RegularTwoPoint: return new RegularCalibration(table);
    case SecondOrderPoly: return new SecondOrderPolyCalibration(table);
    case ThirdOrderPoly: return new ThirdOrderPolyCalibration(table);
    case ForthDegreeThreePoint: return new ForthDegreeThreePointCalibration(table);
    case ExpThreePoint: return new ExpThreePointCalibration(table);
    default: return (Calibration*)0x00;
    }
}

Calibration* Calibration::create
(
    CalibrationType type,
    const MassList &ms,
    const PeakTable &ps
)
{
    switch(type)
    {
    case RegularTwoPoint: return new RegularCalibration(ms, ps);
    case SecondOrderPoly: return new SecondOrderPolyCalibration(ms, ps);
    case ThirdOrderPoly: return new ThirdOrderPolyCalibration(ms, ps);
    case ForthDegreeThreePoint: return new ForthDegreeThreePointCalibration(ms, ps);
    case ExpThreePoint: return new ExpThreePointCalibration(ms, ps);
    default: return (Calibration*)0x00;
    }
}

Calibration::CalibrationType Calibration::strToType(const String &strType)
{
    std::map<String, CalibrationType>::iterator it = s_supportedTypes.find(strType);
    if(it != s_supportedTypes.end())
    {
        return it->second;
    }
    else
    {
        return ErrorType;
    }
}

Calibration::StringList Calibration::supportedTypes()
{
    StringList ret;
    for(const auto& it : s_supportedTypes)
    {
        ret.push_back(it.first);
    }
    return ret;
}

Regression::XYWData Calibration::createRegressionData(const MassTimeIntensity &t) const
{
    XYWData ret;
    for(const CalibTabEntry& p : t)
    {
        ret.emplace(trToRegData(p));
    }
    return ret;
}

Regression::XYWData Calibration::createRegressionData
(
    const MassList &ms,
    const PeakTable &ps
) const
{
    assert(ms.size() == ps.size());
    XYWData ret;
    MassList::const_iterator itm = ms.begin();
    PeakTable::const_iterator itp = ps.begin();
    for(; itm != ms.end(); ++itm, ++itp)
    {
        ret.emplace(trToRegData(*itm, *itp));
    }
    return ret;
}

RegularCalibration::RegularCalibration(const MassTimeIntensity &t)
{
    XYWData reg = createRegressionData(t);
    m_pRegression.reset(new Regression(reg, 1));
}

RegularCalibration::RegularCalibration(const MassList &ms, const PeakTable &ps)
{
    XYWData reg = createRegressionData(ms, ps);
    m_pRegression.reset(new Regression(reg, 1));
}

double RegularCalibration::timeToMass(double time) const
{
    double ret = isValid() ? m_pRegression->yValue(time) : nan("");
    return ret*ret;
}

Calibration::CalibrationType RegularCalibration::type() const
{
    return RegularTwoPoint;
}

bool RegularCalibration::isValid() const
{
    return m_pRegression && m_pRegression->isValid();;
}

double RegularCalibration::t0() const
{
    return - m_pRegression->coefs()[0] / m_pRegression->coefs()[1];
}

Regression::XYWPoint RegularCalibration::trToRegData(const CalibTabEntry &p) const
{
    return {p.second.first, {::sqrt(p.first), p.first * p.second.second}};
}

Regression::XYWPoint RegularCalibration::trToRegData(double mz, const Peak &p) const
{
    return {p.first, {::sqrt(mz), p.second}};
}

SecondOrderPolyCalibration::SecondOrderPolyCalibration(const MassTimeIntensity &t)
{
    XYWData reg = createRegressionData(t);
    m_pRegression.reset(new Regression(reg, 2));
}

SecondOrderPolyCalibration::SecondOrderPolyCalibration
(
    const MassList &ms,
    const PeakTable &ps
)
{
    XYWData reg = createRegressionData(ms, ps);
    m_pRegression.reset(new Regression(reg, 2));
}

double SecondOrderPolyCalibration::timeToMass(double time) const
{
    return isValid() ? m_pRegression->yValue(time) : ::nan("");
}

Calibration::CalibrationType SecondOrderPolyCalibration::type() const
{
    return SecondOrderPoly;
}

bool SecondOrderPolyCalibration::isValid() const
{
    return m_pRegression && m_pRegression->isValid();;
}

Regression::XYWPoint SecondOrderPolyCalibration::trToRegData(const CalibTabEntry &p) const
{
    return {p.second.first, {p.first, p.second.second}};
}

Regression::XYWPoint SecondOrderPolyCalibration::trToRegData(double mz, const Peak &p) const
{
    return {p.first, {mz, p.second}};
}

ThirdOrderPolyCalibration::ThirdOrderPolyCalibration(const MassTimeIntensity &t)
{
    XYWData reg = createRegressionData(t);
    m_pRegression.reset(new Regression(reg, 3));
}

ThirdOrderPolyCalibration::ThirdOrderPolyCalibration
(
    const MassList &ms,
    const PeakTable &ps
)
{
    XYWData reg = createRegressionData(ms, ps);
    m_pRegression.reset(new Regression(reg, 3));
}

double ThirdOrderPolyCalibration::timeToMass(double time) const
{
    return isValid() ? m_pRegression->yValue(time) : nan("");
}

Calibration::CalibrationType ThirdOrderPolyCalibration::type() const
{
    return ThirdOrderPoly;
}

bool ThirdOrderPolyCalibration::isValid() const
{
    return m_pRegression && m_pRegression->isValid();;
}

Regression::XYWPoint ThirdOrderPolyCalibration::trToRegData(const CalibTabEntry &p) const
{
    return {p.second.first, {p.first, p.second.second}};
}

Regression::XYWPoint ThirdOrderPolyCalibration::trToRegData(double mz, const Peak &p) const
{
    return {p.first, {mz, p.second}};
}

ForthDegreeThreePointCalibration::ForthDegreeThreePointCalibration
(
    const MassTimeIntensity &table
)
    :
      m_t0(RegularCalibration(table).t0()),
      m_n(0.0)
{

    for(const CalibTabEntry& p : table)
    {
        m_n = std::max(m_n, p.second.first);
    }

    m_t0 /= m_n;

    double t0Next = 1.0;

    size_t nIterLimit = 1000001;
    while (--nIterLimit != 0
           && ::fabs(t0Next - m_t0) > std::numeric_limits<double>::epsilon())
    {
        calculateRegression(table);
        double df = 0.0, f = 0.0;
        for(const CalibTabEntry& p : table)
        {
            double t = p.second.first / m_n - m_t0;
            double d = (2. * m_pRegression->coefs()[0]
                    + 4. * m_pRegression->coefs()[1] * t * t) * t;
            f += d * (m_pRegression->yValue(t) - p.first);
            df += d*d;
        }
        if(df == 0.0) break;
        t0Next = m_t0;
        m_t0 += f / df;
    }
    calculateRegression(table);
}

ForthDegreeThreePointCalibration::ForthDegreeThreePointCalibration
(
    const MassList &ms,
    const PeakTable &ps
)
    :
      m_t0(RegularCalibration(ms,ps).t0()),
      m_n(0.0)
{
    for(const Peak& p : ps)
    {
        m_n = std::max(m_n, p.first);
    }

    m_t0 /= m_n;

    double t0Next = 1.0;

    size_t nIterLimit = 1000001;
    while (--nIterLimit != 0
           && ::fabs(t0Next - m_t0) > std::numeric_limits<double>::epsilon())
    {
        calculateRegression(ms, ps);
        double df = 0.0, f = 0.0;
        MassList::const_iterator itm = ms.begin();
        PeakTable::const_iterator itp = ps.begin();
        for(; itm != ms.end(); ++itm, ++itp)
        {
            double t = itp->first / m_n - m_t0;
            double d = (2. * m_pRegression->coefs()[0]
                    + 4. * m_pRegression->coefs()[1] * t * t) * t;
            f += d * (m_pRegression->yValue(t) - *itm);
            df += d*d;
        }
        if(df == 0.0) break;
        t0Next = m_t0;
        m_t0 += f / df;
    }
    calculateRegression(ms, ps);
}

double ForthDegreeThreePointCalibration::timeToMass(double time) const
{
    return isValid() ? m_pRegression->yValue(time/m_n - m_t0) : ::nan("");
}

Calibration::CalibrationType ForthDegreeThreePointCalibration::type() const
{
    return ForthDegreeThreePoint;
}

bool ForthDegreeThreePointCalibration::isValid() const
{
    return m_pRegression && m_pRegression->isValid();;
}

Regression::XYWPoint ForthDegreeThreePointCalibration::trToRegData
(
    const CalibTabEntry &p
) const
{
    return {p.second.first/m_n - m_t0, {p.first, p.second.second}};
}

Regression::XYWPoint ForthDegreeThreePointCalibration::trToRegData
(
    double mz,
    const Peak &p
) const
{
    return {p.first/m_n - m_t0, {mz, p.second}};
}

LinearRegression::VectorFun ForthDegreeThreePointCalibration::basisFuns()
{
    return VectorFun
    {
        [](double x)->double {return x*x;},
        [](double x)->double
        {
            double x2 = x*x;
            return x2*x2;
        }
    };
}

void ForthDegreeThreePointCalibration::calculateRegression(const MassTimeIntensity &t)
{
    XYWData reg = createRegressionData(t);
    m_pRegression.reset(new LinearRegression(reg, basisFuns()));
}

void ForthDegreeThreePointCalibration::calculateRegression
(
    const MassList &ms,
    const PeakTable &ps
)
{
    XYWData reg = createRegressionData(ms, ps);
    m_pRegression.reset(new LinearRegression(reg, basisFuns()));
}

ExpThreePointCalibration::ExpThreePointCalibration(const MassTimeIntensity &table)
    :
    m_t0(RegularCalibration(table).t0())
{
    double n = 0.0;
    for(const CalibTabEntry& p : table)
    {
        n = std::max(n, p.second.first);
    }

    double t0Next = n;

    size_t nIterLimit = 1000001;
    while (--nIterLimit != 0
           && ::fabs((t0Next - m_t0)/n) > std::numeric_limits<double>::epsilon())
    {
        calculateRegression(table);
        double df = 0.0, f = 0.0;
        for(const CalibTabEntry& p : table)
        {
            double t = p.second.first - m_t0;
            double d = m_pRegression->coefs()[1] / t;
            f += d * (m_pRegression->yValue(::log(t)) - ::log(p.first));
            df += d*d;
        }
        if(df == 0.0) break;
        t0Next = m_t0;
        m_t0 += f / df;
    }
    calculateRegression(table);
}

ExpThreePointCalibration::ExpThreePointCalibration(const MassList &ms, const PeakTable &ps)
    :
      m_t0(RegularCalibration(ms, ps).t0())
{
    double n = 0.0;
    for(const Peak& p : ps)
    {
        n = std::max(n, p.first);
    }

    double t0Next = n;

    size_t nIterLimit = 1000001;
    while (--nIterLimit != 0
           && ::fabs((t0Next - m_t0)/n) > std::numeric_limits<double>::epsilon())
    {
        calculateRegression(ms, ps);
        double df = 0.0, f = 0.0;
        MassList::const_iterator itm = ms.begin();
        PeakTable::const_iterator itp = ps.begin();
        for(; itm != ms.end(); ++itm, ++itp)
        {
            double t = itp->first - m_t0;
            double d = m_pRegression->coefs()[1] / t;
            f += d * (m_pRegression->yValue(::log(t)) - ::log(*itm));
            df += d*d;
        }
        if(df == 0.0) break;
        t0Next = m_t0;
        m_t0 += f / df;
    }
    calculateRegression(ms, ps);
}

double ExpThreePointCalibration::timeToMass(double time) const
{
    return isValid() ? ::exp(m_pRegression->yValue(::log(time - m_t0))) : nan("");
}

Calibration::CalibrationType ExpThreePointCalibration::type() const
{
    return ExpThreePoint;
}

bool ExpThreePointCalibration::isValid() const
{
    return m_pRegression && m_pRegression->isValid();
}

Regression::XYWPoint ExpThreePointCalibration::trToRegData(const CalibTabEntry &p) const
{
    double t = p.second.first - m_t0;
    assert(t > 0);
    return {::log(t), {::log(p.first), p.first * p.second.second}};
}

Regression::XYWPoint ExpThreePointCalibration::trToRegData(double mz, const Peak &p) const
{
    double t = p.first - m_t0;
    assert(t > 0);
    return {::log(t), {::log(mz), mz * p.second}};
}

void ExpThreePointCalibration::calculateRegression(const MassTimeIntensity &t)
{
    XYWData reg = createRegressionData(t);
    m_pRegression.reset(new Regression(reg, 1));
}

void ExpThreePointCalibration::calculateRegression(const MassList &ms, const PeakTable &ps)
{
    XYWData reg = createRegressionData(ms, ps);
    m_pRegression.reset(new Regression(reg, 1));
}
