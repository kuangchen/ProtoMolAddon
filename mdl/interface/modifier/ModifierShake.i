%module ModifierShake
%{
#include <protomol/modifier/Modifier.h>
#include <protomol/modifier/ModifierShake.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
using namespace ProtoMol;
%}

%include <protomol/type/Real.h>
%include <protomol/modifier/Modifier.h>
%include <protomol/modifier/ModifierMetaRattleShake.h>
%include <protomol/modifier/ModifierMetaShake.h>
%include <protomol/modifier/ModifierShake.h>
