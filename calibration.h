#ifndef CALIBRATION_H
#define CALIBRATION_H

#include <set>
#include <map>
#include <string>
#include <memory>
#include <vector>

#include "regression.h"

class Calibration
{
public:
    /**
     * @brief Peak in a time scale {time, Intensity}
     */
    typedef Regression::XYWPoint XYWPoint;
    typedef Regression::XYWData XYWData;
    typedef std::set<double> MassList;
    typedef std::pair<double, double> Peak;
    typedef std::map<double, double> PeakTable;
    typedef std::pair<double, Peak> CalibTabEntry;
    typedef std::map<double, Peak> MassTimeIntensity;
    typedef std::string String;
    typedef std::vector<String> StringList;

    enum CalibrationType
    {
        RegularTwoPoint,
        SecondOrderPoly,
        ThirdOrderPoly,
        ForthDegreeThreePoint,
        ExpThreePoint,
        ErrorType
    };

    typedef std::pair<String, CalibrationType> CalTypeEntry;
    typedef std::map<String, CalibrationType> CalTypesList;

    static Calibration* create(CalibrationType type, const MassTimeIntensity& table);

    static Calibration* create
    (
        CalibrationType type,
        const MassList &ms,
        const PeakTable &ps
    );

    static CalibrationType strToType(const String& strType);

    static StringList supportedTypes();

    virtual CalibrationType type() const = 0;

    virtual bool isValid() const = 0;

    virtual ~Calibration(){}

    /**
     * @brief timeToMass
     * @param time
     * @return Mass value correspondent to a given time value
     */
    virtual double timeToMass(double time) const = 0;

protected:

    /**
     * @brief trToRegData transforms mass peak to a regression point
     * @return regression point
     */
    virtual XYWPoint trToRegData(const CalibTabEntry& p) const = 0;
    virtual XYWPoint trToRegData(double mz, const Peak& p) const = 0;

    XYWData createRegressionData(const MassTimeIntensity& t) const;

    XYWData createRegressionData(const MassList& ms, const PeakTable& ps) const;

private:
    static CalTypesList s_supportedTypes;
};

/**
 * @brief The RegularCalibration class implements regular TOF calibration
 * using equation m = a*(t-t0)*(t-t0);
 */
class RegularCalibration : public Calibration
{
    typedef std::unique_ptr<Regression> PRegression;
public:
    RegularCalibration(const MassTimeIntensity& t);

    RegularCalibration(const MassList& ms, const PeakTable& ps);

    virtual double timeToMass(double time) const;

    virtual CalibrationType type() const;

    virtual bool isValid() const;

    /**
     * @brief t0
     * @return time shift parameter of calibration
     */
    double t0() const;

protected:

    virtual XYWPoint trToRegData(const CalibTabEntry &p) const;

    virtual XYWPoint trToRegData(double mz, const Peak& p) const;

private:
    PRegression m_pRegression;
};

/**
 * @brief The SecondOrderPolyCalibration class implements TOF calibration
 * using equation m = at^2 + bt + c
 */
class SecondOrderPolyCalibration : public Calibration
{
    typedef std::unique_ptr<Regression> PRegression;
public:
    SecondOrderPolyCalibration(const MassTimeIntensity& t);

    SecondOrderPolyCalibration(const MassList& ms, const PeakTable& ps);

    virtual double timeToMass(double time) const;

    virtual CalibrationType type() const;

    virtual bool isValid() const;

protected:

    virtual XYWPoint trToRegData(const CalibTabEntry &p) const;
    virtual XYWPoint trToRegData(double mz, const Peak& p) const;

private:
    PRegression m_pRegression;
};

/**
 * @brief The ThirdOrderPolyCalibration class implements TOF calibration
 * using equation m = at^3 + bt^2 + ct + d
 */
class ThirdOrderPolyCalibration : public Calibration
{
    typedef std::unique_ptr<Regression> PRegression;

public:
    ThirdOrderPolyCalibration(const MassTimeIntensity& t);

    ThirdOrderPolyCalibration(const MassList& ms, const PeakTable& ps);

    virtual double timeToMass(double time) const;

    virtual CalibrationType type() const;

    virtual bool isValid() const;

protected:

    virtual XYWPoint trToRegData(const CalibTabEntry &p) const;

    virtual XYWPoint trToRegData(double mz, const Peak& p) const;

private:
    PRegression m_pRegression;


};

/**
 * @brief The ForthDegreeThreePointCalibrattion class implements TOF calibration
 * using equation m = k2(t-t0)^2 + k4(t-t0)^4
 */
class ForthDegreeThreePointCalibration : public Calibration
{
    typedef std::unique_ptr<LinearRegression> PRegression;
    typedef LinearRegression::Vector Vector;
    typedef LinearRegression::VectorFun VectorFun;
public:

    ForthDegreeThreePointCalibration(const MassTimeIntensity& table);

    ForthDegreeThreePointCalibration(const MassList& ms, const PeakTable& ps);

    virtual double timeToMass(double time) const;

    virtual CalibrationType type() const;

    virtual bool isValid() const;

protected:

    virtual XYWPoint trToRegData(const CalibTabEntry &p) const;

    virtual XYWPoint trToRegData(double mz, const Peak& p) const;

private:

    double m_t0;

    /**
     * @brief m_n notmalizes t values to [0,1] interval
     */
    double m_n;

    PRegression m_pRegression;

    /**
     * @brief basisFuns creates a vector of regression basis functions
     * @return
     */
    static VectorFun basisFuns();

    void calculateRegression(const MassTimeIntensity& t);
    void calculateRegression(const MassList& ms, const PeakTable& ps);
};

/**
 * @brief The ExpThreePointCalibration class implements TOF calibration
 * using equation m = kb (t-t0)^b
 */
class ExpThreePointCalibration : public Calibration
{
    typedef std::unique_ptr<Regression> PRegression;
    typedef Regression::Vector Vector;
public:
    ExpThreePointCalibration(const MassTimeIntensity& table);

    ExpThreePointCalibration(const MassList& ms, const PeakTable& ps);

    virtual double timeToMass(double time) const;

    virtual CalibrationType type() const;

    virtual bool isValid() const;

protected:

    virtual XYWPoint trToRegData(const CalibTabEntry &p) const;
    virtual XYWPoint trToRegData(double mz, const Peak& p) const;

private:
    double m_t0;

    PRegression m_pRegression;

    void calculateRegression(const MassTimeIntensity& t);
    void calculateRegression(const MassList& ms, const PeakTable& ps);
};

#endif // CALIBRATION_H
