#pragma once

#include <concepts>

template <typename T>
concept SupportClosedUnaryMinus = requires(const T& val) {
  { -val } -> std::same_as<T>;
};

template <typename T, typename Q>
concept AddableWith = requires(const T& lhv, const Q& rhv) { lhv + rhv; };

template <typename T>
concept Addable = AddableWith<T, T>;

template <typename T>
concept AddableClosed = Addable<T> and requires(const T& lhv, const T& rhv) {
  { lhv + rhv } -> std::same_as<T>;
};

template <typename T, typename Q>
concept SubstractableWith = requires(const T& lhv, const Q& rhv) { lhv - rhv; };

template <typename T>
concept Substractable = SubstractableWith<T, T>;

template <typename T>
concept SubstractableClosed =
    Substractable<T> and requires(const T& lhv, const T& rhv) {
      { lhv - rhv } -> std::same_as<T>;
    };

template <typename T, typename Q>
concept MultipliableWith = requires(const T& lhv, const Q& rhv) { lhv * rhv; };

template <typename T>
concept Multipliable = MultipliableWith<T, T>;

template <typename T>
concept MultipliableClosed =
    Multipliable<T> and requires(const T& lhv, const T& rhv) {
      { lhv * rhv } -> std::same_as<T>;
    };

template <typename T, typename Q>
concept DivisibleWith = requires(const T& lhv, const Q& rhv) { lhv / rhv; };

template <typename T>
concept Divisible = DivisibleWith<T, T>;

template <typename T>
concept DivisibleClosed =
    Divisible<T> and requires(const T& lhv, const T& rhv) {
      { lhv / rhv } -> std::same_as<T>;
    };

template <typename T>
concept Arithemtical = AddableClosed<T> and SubstractableClosed<T> and
                       DivisibleClosed<T> and MultipliableClosed<T>;
