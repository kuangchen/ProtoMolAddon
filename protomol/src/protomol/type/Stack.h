/* -*- c++ -*- */
#ifndef STACK_H
#define STACK_H


#include <string.h>

namespace ProtoMol {
  //________________________________________ Stack
  /**
   * Stack which also enables random access.
   */
  template<class T>
  class Stack {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    Stack();
    ~Stack();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Stack
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    void addElement(T newElement);
    T popElement();
    T getElement(unsigned int index);
    unsigned int getNumElements();
    void reset(bool delmem);
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    T *m_elem;
    unsigned int m_numElements;
    unsigned int m_size;
    static const unsigned int INCR_SIZE = 50;
  };

  //________________________________________ INLINES
  template<class T>
  Stack<T>::Stack() {
    m_elem = NULL;
    m_numElements = 0;
    m_size = 0;
  }

  template<class T>
  Stack<T>::~Stack() {
    if (m_elem != NULL)
      delete[] m_elem;
    m_elem = NULL;
  }

  template<class T>
  void Stack<T>::addElement(T newElement) {
    T *temp;

    if (m_numElements >= m_size) {
      temp = new T[m_size + INCR_SIZE];
      if (m_numElements > 0) {
        memcpy(temp, m_elem, sizeof(T) * m_numElements);
        delete[] m_elem;
      }
      m_size += INCR_SIZE;
      m_elem = temp;
    }
    m_elem[m_numElements] = newElement;
    m_numElements++;
  }

  template<class T>
  T Stack<T>::popElement() {
    T result;

    result = m_elem[m_numElements - 1];
    m_numElements--;

    return result;
  }

  template<class T>
  T Stack<T>::getElement(unsigned int index) {
    T result = T();

    if (index < m_numElements)
      result = m_elem[index];

    return result;
  }

  template<class T>
  void Stack<T>::reset(bool delMem) {
    if (delMem == true) {
      if (m_size > 0)
        delete[] m_elem;
      m_elem = NULL;
      m_size = 0;
      m_numElements = 0;
    } else
      m_numElements = 0;
  }

  template<class T>
  unsigned int Stack<T>::getNumElements() {
    return m_numElements;
  }
}

#endif
