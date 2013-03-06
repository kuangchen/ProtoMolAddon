#ifndef INPUT_VALUE_MACROS_H
#define INPUT_VALUE_MACROS_H

#include <string>
#include <vector>

#include <protomol/type/Vector.h>
#include <protomol/config/ConstraintValueType.h>
#include <protomol/config/ValueType.h>

#define declareInputValue(NAME, TYPE, CONSTRAINT)                        \
  struct NAME##Identifier {                                              \
    virtual ~NAME##Identifier() {}                                       \
    static const std::string keyword;                                    \
    static const std::vector<std::string> aliases;                       \
    static const std::string text;                                       \
  };                                                                     \
  template <typename B, ValueType::Enum T, ConstraintValueType::Enum C>  \
  class InputValue;                                                      \
  typedef InputValue<NAME##Identifier, ValueType::TYPE,                  \
                     ConstraintValueType::CONSTRAINT> NAME;

#define defineInputValue(NAME, KEYWORD)                                  \
  const string NAME##Identifier::keyword (KEYWORD);                      \
  const vector<string> NAME##Identifier::aliases;                        \
  const string NAME##Identifier::text("");

#define defineInputValueWithAliases(NAME, KEYWORD, ALIASES)              \
  const string NAME##Identifier::keyword(KEYWORD);                       \
  const vector<string> NAME##Identifier::                                \
    aliases(static_cast<vector<string> >(Vector<string>ALIASES));        \
  const string NAME##Identifier::text("");

#define defineInputValueAndText(NAME, KEYWORD, TEXT)                     \
  const string NAME##Identifier::keyword(KEYWORD);                       \
  const vector<string> NAME##Identifier::aliases;                        \
  const string NAME##Identifier::text(TEXT);

#define defineInputValueWithAliasesAndText(NAME, KEYWORD, ALIASES, TEXT) \
  const string NAME##Identifier::keyword(KEYWORD);                       \
  const vector<string> NAME##Identifier::                                \
    aliases(static_cast<vector<string> >(Vector<string>ALIASES));        \
  const string NAME##Identifier::text(TEXT);

#endif // INPUT_VALUE_MACROS_H
