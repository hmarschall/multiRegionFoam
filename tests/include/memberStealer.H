/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    A set of macro to ease testing of private members/methods of classes

Author
    Mohammed Elwardi Fadeli (elwardifadeli@gmail.com)

Notes
    Not recommended to use with your own classes; make the methods public
    instead
    Destined to test behaviour of FE/OF classes

\*---------------------------------------------------------------------------*/

#ifndef memberStealer_H
#define memberStealer_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Example Usage:

// In global space / outer scope:

    #include "memberStealer.H"

    // Example Class
    class MyClass
    {
    public:
        MyClass(){};
        virtual ~MyClass(){};
    
        int getValue() { return value; }
    
    private:
        int value = 0;
        int setValue(int a) { value = a; }
    };
    
    // Explicit specialization for the private member in example class
    SPECIALIZE_MEMBER_STEALER(value, int, MyClass);

    // Explicit specialization for the private member method
    using MethodType = int(MyClass::*)(int);
    SPECIALIZE_MEMBER_METHOD_STEALER(setValue, MethodType, MyClass);


// Where you want to use the private member or member method

    MyClass obj; // Example object of example type
    // Get a ref to member variable
    auto& number = GET_REF(value, MyClass, obj);
    // Call member method setValue with args (10)
    CALL_MEMBER_METHOD(setValue, MyClass, obj) (10);

 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// A template to generate steal structs
template<typename Tag, typename Tag::type M>
struct memberStealer
{
    friend typename Tag::type get(Tag) {
        return M;
    }
};

// Tag used to forcibly gain access to memberClass::memberName
// This is an explicit specialization
#define SPECIALIZE_MEMBER_STEALER(memberName, memberType, memberClass)        \
struct memberName##TagIN##memberClass                                         \
{                                                                             \
    typedef memberType memberClass::*type;                                    \
    friend type get(memberName##TagIN##memberClass);                          \
};                                                                            \
template struct memberStealer                                                 \
<                                                                             \
    memberName##TagIN##memberClass,                                           \
    &memberClass::memberName                                                  \
>

// Tag used to forcibly gain access to memberClass::*
// This is an explicit specialization
#define SPECIALIZE_MEMBER_METHOD_STEALER(memberName, memberType, memberClass) \
struct memberName##TagIN##memberClass                                         \
{                                                                             \
    typedef memberType type;                                                  \
    friend type get(memberName##TagIN##memberClass);                          \
};                                                                            \
template struct memberStealer                                                 \
<                                                                             \
    memberName##TagIN##memberClass,                                           \
    &memberClass::memberName                                                  \
>

#define GET_REF(memberName, memberClass, object)                              \
object.*get(memberName##TagIN##memberClass())

#define CALL_MEMBER_METHOD(memberName, memberClass, object)                   \
(object.*get(memberName##TagIN##memberClass()))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
