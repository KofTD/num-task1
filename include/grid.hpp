#pragma once

#include <concepts>
#include <cstdlib>
#include <format>
#include <limits>
#include <stdexcept>
#include <utility>

#include "fp-comprasion.hpp"

template <std::floating_point F>
class UniformGrid {
 private:
  F start_;
  F end_;
  F step_;

 public:
  UniformGrid() = delete;
  UniformGrid(const UniformGrid&) = default;
  UniformGrid(UniformGrid&&) = default;
  auto operator=(const UniformGrid&) -> UniformGrid& = default;
  auto operator=(UniformGrid&&) -> UniformGrid& = default;

  UniformGrid(F start, F end, F step) : start_(start), end_(end), step_(step) {
    if (start > end) {
      throw std::invalid_argument("Start must be less than end");
    }
  };
  UniformGrid(F start, F end, std::unsigned_integral auto steps)
      : start_(start),
        end_(end),
        step_(static_cast<F>((end_ - start_) / steps)) {
    if (start > end) {
      throw std::invalid_argument("Start must be less than end");
    }
  }

  class Iterator {
   private:
    const UniformGrid* grid_ = nullptr;
    std::size_t step_ = 0;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = std::pair<F, F>;
    using reference = value_type;

    Iterator() = default;
    Iterator(const Iterator&) = default;
    Iterator(Iterator&&) = default;
    auto operator=(const Iterator&) -> Iterator& = default;
    auto operator=(Iterator&&) -> Iterator& = default;
    explicit Iterator(UniformGrid* grid, std::size_t step = 0)
        : grid_(grid), step_(step) {}

    auto operator*() const -> value_type { return grid_->region(step_); }
    auto operator++() -> Iterator& {
      step_++;
      return *this;
    }
    auto operator++(int) -> Iterator {
      auto tmp = *this;
      ++*this;
      return tmp;
    }
    auto operator--() -> Iterator& {
      step_--;
      return *this;
    }
    auto operator--(int) -> Iterator {
      auto tmp = *this;
      --*this;
      return tmp;
    }
    auto operator+=(difference_type diff) -> Iterator& {
      step_ += diff;
      return *this;
    }
    auto operator-=(difference_type diff) -> Iterator& {
      step_ -= diff;
      return *this;
    }
    friend auto operator+(Iterator it, difference_type diff) -> Iterator {
      it += diff;
      return it;
    }
    friend auto operator-(Iterator it, difference_type diff) -> Iterator {
      it -= diff;
      return it;
    }
    friend auto operator-(const Iterator& lhs, const Iterator& rhs)
        -> difference_type {
      return static_cast<difference_type>(rhs.step_) -
             static_cast<difference_type>(lhs.step_);
    }
    auto operator==(const Iterator& other) const -> bool = default;
    auto operator!=(const Iterator& other) const -> bool = default;
    auto operator<=>(const Iterator& other) const {
      return step_ <=> other.step_;
    }
    ~Iterator() = default;
  };
  auto region(std::size_t index) const -> std::pair<F, F> {
    if (index >= steps()) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps()));
    }
    return {node(index), node(index + 1)};
  }
  auto node(std::size_t index) const -> F {
    if (index >= steps()) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps()));
    }
    return start_ + (step_ * index);
  }
  auto step() const noexcept -> F { return step_; }
  [[nodiscard]] auto steps() const noexcept -> std::size_t {
    return static_cast<std::size_t>((end_ - start_) / step_);
  }
  auto begin() -> Iterator { return Iterator(this); }
  auto end() -> Iterator { return Iterator(this, steps()); }
  auto begin() const -> Iterator { return {this}; }
  auto end() const -> Iterator { return {this, steps()}; }

  ~UniformGrid() = default;
};

template <std::floating_point F>
auto midpoint(std::pair<F, F> bounds) -> F {
  return (bounds.second - bounds.first) / 2;
}
