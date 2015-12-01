#ifndef __SBB_MINMAX_HH__
#define __SBB_MINMAX_HH__

#include <vector>
#include "randomvariable.hh"

namespace sbb {

/** Implements the density of a random variable being the maximum of several independent random
 * variables.
 * @ingroup internal */
class MaximumDensityObj: public DensityObj
{
public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MaximumDensityObj(const std::vector<VarObj *> &variables);
  /** Destructor. */
  virtual ~MaximumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  /** Returns the vector of densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
};


/** Implements the density of a random variable being the minimum of several independent random
 * variables.
 * @ingroup internal */
class MinimumDensityObj: public DensityObj
{
public:
  /** Constructor.
   * @param variables Specifies the vector of independent random variables.
   * @throws AssumtionError If the given random variables are not mutually independent. */
  MinimumDensityObj(const std::vector<VarObj *> &variables);
  /** Destructor. */
  virtual ~MinimumDensityObj();
  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void sample(Eigen::VectorXd &out) const;

  /** Returns the vector of densities. */
  inline const std::vector<DensityObj *> &densities() const { return _densities; }

protected:
  /** The vector of densities. */
  std::vector<DensityObj *> _densities;
};


/** Implements the random variable being the maximum of several independent random
 * variables.
 * @ingroup internal */
class MaximumObj: public VarObj
{
public:
  /** Constructor from two random variables.
   * @throws AssumptionError If these two random variables are not independent. */
  MaximumObj(VarObj *a, VarObj *b, const std::string &name="");
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MaximumObj(const std::vector<VarObj *> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MaximumObj();
  virtual void mark();

  virtual DensityObj *density();

  /** Returns the vector of random variables, this RV depends on. */
  const std::vector<VarObj *> variables() const {
    return _variables;
  }

protected:
  /** The vector of random variables, this RV depends on. */
  std::vector<VarObj *> _variables;
  /** The density object. */
  MaximumDensityObj *_density;
};


/** Implements the random variable being the minimum of several independent random
 * variables.
 * @ingroup internal */
class MinimumObj: public VarObj
{
public:
  /** Constructor from two random variables.
   * @throws AssumptionError If these two random variables are not independent. */
  MinimumObj(VarObj *a, VarObj *b, const std::string &name="");
  /** Constructor from a vector of random variables.
   * @throws AssumptionError If these random variables are not independent. */
  MinimumObj(const std::vector<VarObj *> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~MinimumObj();

  virtual void mark();

  virtual DensityObj *density();

  /** Returns the vector of random variables, this RV depends on. */
  const std::vector<VarObj *> variables() const {
    return _variables;
  }

protected:
  /** The vector of random variables, this RV depends on. */
  std::vector<VarObj *> _variables;
  /** The density object of this random variable. */
  MinimumDensityObj *_density;
};

}

#endif // __SBB_MINMAX_HH__
