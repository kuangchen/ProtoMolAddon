/*  -*- c++ -*-  */
#ifndef INPUTVALUE_H
#define INPUTVALUE_H

#include <protomol/config/Value.h>
#include <protomol/config/InputValueMacros.h>
#include <protomol/config/Configuration.h>

#include <protomol/base/Report.h>

namespace ProtoMol {
  //________________________________________________________ InputValue
  /**
     Template class with macros for simple creation of input parameters
     for the configuration file.@n @n

     Example for the temperatur, real number, non-negative @n @n

     // Declaration @n
     declareInputValue(InputTemperature,REAL,NOTNEGATIVE); @n @n

     // Definition, with and without aliases @n
     defineInputValue(InputTemperature,"temperature"); @n
     defineInputValueWithAliases(InputTemperature,"temperature",
       ("temp")("t")); @n @n


     // Defintion with help text @n
     defineInputValueAndText(InputTemperature,"temperature",
       "initial temperature, [K]"); @n
     defineInputValueWithAliasesAndText(InputTemperature,"temperature",
       ("temp")("t"),"initial temperature, [K]"); @n @n

     // Usage @n
     InputTemperature t(9.9); @n
     t = -3;    \\ now t invalid  @n
     t = "3.3"; \\ ok again @n
   */
  template<typename TBase, ValueType::Enum type,
           ConstraintValueType::Enum constraint =
           ConstraintValueType::NOCONSTRAINTS>
  class InputValue : public TBase {
  public:
    typedef typename Enum2ValueTraits<type>::Type Type;
    typedef typename ConstraintValueEnum::Enum2Type<constraint> Constraint;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    InputValue() : myValue(typename Type::Type(), Constraint(),
                           Value::undefined) {}
    template<typename T>
    explicit InputValue(const T &v) : myValue(Type(), Constraint(), v, true) {}
    explicit InputValue(const char *v) :
      myValue(Type(), Constraint(), std::string(v), true) {}
    virtual ~InputValue() {}

    // Assignment
    template<typename T>
    InputValue &operator=(const T &value) {
      myValue.set(value);
      return *this;
    }

    InputValue &operator=(const char *value) {
      return operator=(std::string (value));
    }

    InputValue &operator=(const InputValue &rhs) {
      InputValue(rhs).swap(*this);
      return *this;
    }

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class InputValue
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    operator Value() const {return myValue;}

    template<typename T>
    operator T() const {
      T val = T();
      myValue.get(val);
      return val;
    }

    template<typename T>
    bool get(T &val) const {return myValue.get(val);}

    bool valid() const {return myValue.valid();}

    InputValue &swap(InputValue &rhs) {
      std::swap(myValue, rhs.myValue);
      return *this;
    }

    static void registerConfiguration(Configuration *config, Value v) {
      config->registerKeyword(TBase::keyword, v);
      config->registerAliases(TBase::keyword, TBase::aliases);
      if (!TBase::text.empty())
        config->setText(TBase::keyword, TBase::text);
    }

    template<typename T>
    static void registerConfiguration(Configuration *config, T v) {
      config->registerKeyword(TBase::keyword, Value(v));
      config->registerAliases(TBase::keyword, TBase::aliases);
      if (!TBase::text.empty())
        config->setText(TBase::keyword, TBase::text);
    }

    static void registerConfiguration(Configuration *config) {
      config->registerKeyword(TBase::keyword,
        Value(Type::empty(), Constraint(), Value::undefined));
      config->registerAliases(TBase::keyword, TBase::aliases);
      if (!TBase::text.empty())
        config->setText(TBase::keyword, TBase::text);
    }

    static void registerConfiguration(Configuration *config,
                                      Value v,
                                      const Text &txt) {
      config->registerKeyword(TBase::keyword, v);
      config->registerAliases(TBase::keyword, TBase::aliases);
      config->setText(TBase::keyword, txt.text);
    }

    template<typename T>
    static void registerConfiguration(Configuration *config,
                                      T v,
                                      const Text &txt) {
      config->registerKeyword(TBase::keyword, Value(v));
      config->registerAliases(TBase::keyword, TBase::aliases);
      config->setText(TBase::keyword, txt.text);
    }

    static void registerConfiguration(Configuration *config, const Text &txt) {
      config->registerKeyword(TBase::keyword,
        Value(Type::empty(), Constraint(), Value::undefined));
      config->registerAliases(TBase::keyword, TBase::aliases);
      config->setText(TBase::keyword, txt.text);
    }

    friend Report::MyStreamer &operator<<(Report::MyStreamer &OS,
                                          const InputValue &v) {
      OS << v.myValue.getString();
      return OS;
    }


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // private data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Value myValue;
  };

};
#endif
