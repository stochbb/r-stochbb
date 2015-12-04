#ifndef CHAIN_HH
#define CHAIN_HH

#include <vector>
#include "api.hh"
#include "randomvariable.hh"

namespace sbb {

/** Implements the convolution of several PDFs using the FFT convolution. */
class ConvolutionDensityObj: public DensityObj
{
protected:
  /** Constructs a new PDF as the convolution of the PDFs of the given variables
   * (weak references). */
  ConvolutionDensityObj(const std::vector<VarObj *> &variables);

public:
  /** Constructs a new PDF as the convolution of the PDFs of the given variables. */
  ConvolutionDensityObj(const std::vector<Var> &variables);
  /** Destructor. */
  virtual ~ConvolutionDensityObj();

  virtual void mark();

  virtual void eval(double Tmin, double Tmax, Eigen::VectorXd &out) const;
  virtual void evalCDF(double Tmin, double Tmax, Eigen::VectorXd &out) const;

  /** Returns the number of underlaying densities. */
  inline size_t numDensities() const { return _densities.size(); }

  /** Returns the i-th density. */
  inline Density density(size_t i) const {
    _densities[i]->ref();
    return _densities[i];
  }

  /** Returns a vector of weak references to the densities of the underlaying variables. */
  inline const std::vector<DensityObj *> &densities() const {
    return _densities;
  }

  int compare(const DensityObj &other) const;
  void print(std::ostream &stream) const;

protected:
  /** Inernal used function to perfrom some convolutions analytically. */
  void _combine_densities();

protected:
  /** The @c DensityObj instances of the underlaying variables. */
  std::vector<DensityObj *> _densities;
  friend class ChainObj;
};


/** Represetns the sum of several independent random variables. */
class ChainObj : public DerivedVarObj
{
public:
  /** Constructs the sum of the given random variables. */
  ChainObj(const std::vector<Var> &variables, const std::string &name="");
  /** Destructor. */
  virtual ~ChainObj();
  virtual void mark();

  virtual Density density();


protected:
  /** The density of the sum, the convolution of all PDFs of the underlaying random variables. */
  ConvolutionDensityObj *_density;
};

}

#endif // CHAIN_HH
