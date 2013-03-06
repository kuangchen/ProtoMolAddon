#ifndef SINGLETON_H
#define SINGLETON_H

template< typename Type >
class Singleton
{
    private:
        static Type* pInstance;

    public:
        static Type& Instance() {
            if ( Singleton<Type>::pInstance == 0 ) {
                Singleton<Type>::pInstance = new Type;
            }

            return *Singleton<Type>::pInstance;
        }
};

template< typename Type > Type* Singleton<Type>::pInstance = 0;

#endif // SINGLETON_H
