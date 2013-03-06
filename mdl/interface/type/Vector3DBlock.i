%module Vector3DBlock
%{
#include <protomol/type/Vector3DBlock.h>
#include <protomol/type/Real.h>
#include <protomol/ProtoMolApp.h>
#include "ndarrayobject.h"
using namespace ProtoMol;
%}

%include <protomol/type/Vector3DBlock.h>

%extend ProtoMol::Vector3DBlock {


   void delC() {
      if (self->c) delete self->c;
   }

   Real __getitem__(int index){
	return self->c[index];
   }

   void __setitem__(int index, Real val) {
	self->c[index] = val;
   }


   void printC() {
      cout << "C: " << self->c << " MEMBER C: " << (*self)[0].c << endl;
   }

   void setC(PyObject* rhs) {
      import_array();
      //self->resize(((PyArrayObject*)rhs)->dimensions[0]/3);
      //if ((PyArrayObject*)(self->c)) {
      if (self->c != NULL && (self->c != (Real*)(((PyArrayObject*)rhs)->data))) {
	delete self->c;
        self->c = 0;
        self->c = (Real*)(((PyArrayObject*)rhs)->data);
        int n = ((PyArrayObject*)rhs)->dimensions[0]/3;
      if (n < self->size()) {
         for (unsigned int i = 0; i < n; i++) {
                self->vec[i].c = self->c+3*i;
         }
         self->vec.resize(n);

      }
      else if (n > self->size()) {
         for (unsigned int i = 0; i < self->size(); i++) {
		self->vec[i].c = self->c+3*i;
	 }								        
         for (unsigned int i = self->size(); i < n; i++) {
            self->vec.push_back(Vector3DB(self->c[3*i], self->c[3*i+1], self->c[3*i+2], (self->c)+3*i));
	 }
      }
      else {
        for (unsigned int i = 0; i < n; i++) {
	       self->vec[i].c = self->c+3*i;
	}
      }
      //(*self)[0].c = self->c;
      }
   }

   PyObject* getC() {
	int mySize = self->size()*3;
        npy_intp dims[1] = {mySize};
	import_array1(NULL);
	PyObject* rhs = PyArray_SimpleNewFromData(1,dims,PyArray_DOUBLE,(char*)(self->c));
        return rhs;
   }

};
