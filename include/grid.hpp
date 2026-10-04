#pragma once

#include <concepts>
#include <cstdlib>
#include <format>
#include <stdexcept>
#include <utility>

template <std::floating_point F>
class UniformGrid {
 private:
  F start_;
  F end_;
  F step_;
  std::size_t steps_;

 public:
  UniformGrid() = delete;
  UniformGrid(const UniformGrid&) = default;
  UniformGrid(UniformGrid&&) = default;
  auto operator=(const UniformGrid&) -> UniformGrid& = default;
  auto operator=(UniformGrid&&) -> UniformGrid& = default;

  template <typename S>
    requires std::convertible_to<S, std::size_t>
  UniformGrid(F start, F end, S steps)
      : start_(start), end_(end), steps_(steps), step_((end - start) / steps) {
    if (start > end) {
      throw std::invalid_argument("Start must be less than end");
    }
    if (steps <= 0) {
      throw std::invalid_argument("Steps must be a positive number");
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
    friend auto operator+(Iterator iter, difference_type diff) -> Iterator {
      iter += diff;
      return iter;
    }
    friend auto operator-(Iterator iter, difference_type diff) -> Iterator {
      iter -= diff;
      return iter;
    }
    friend auto operator-(const Iterator& lhs, const Iterator& rhs)
        -> difference_type {
      return static_cast<difference_type>(rhs.step_) -
             static_cast<difference_type>(lhs.step_);
    }
    auto operator==(const Iterator& other) const -> bool = default;
    auto operator<=>(const Iterator& other) const {
      return step_ <=> other.step_;
    }
    ~Iterator() = default;
  };
  auto region(std::size_t index) const -> std::pair<F, F> {
    if (index > steps_) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps_));
    }
    return {node(index), node(index + 1)};
  }
  auto node(std::size_t index) const -> F {
    if (index > steps_) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps_));
    }
    return start_ + (step_ * index);
  }

  auto supNode(std::size_t index) const -> F {
    if (index > steps_ - 1) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps_ - 1));
    }
    auto main_node = this->node(index);
    return main_node + (step_ / 2);
  }

  auto supRegion(std::size_t index) const -> std::pair<F, F> {
    if (index > steps_ - 2) {
      throw std::invalid_argument(
          std::format("Index {} is out of bounds [0, {}]", index, steps_ - 2));
    }
    auto left = this->supNode(index);
    auto right = this->supNode(index + 1);
    return {left, right};
  }

  auto step() const noexcept -> F { return step_; }
  [[nodiscard]] auto steps() const noexcept -> std::size_t { return steps_; }
  [[nodiscard]] auto nodes() const noexcept -> std::size_t {
    return steps_ + 1;
  }
  auto begin() -> Iterator { return Iterator(this); }
  auto end() -> Iterator { return Iterator(this, steps()); }
  auto begin() const -> Iterator { return {this}; }
  auto end() const -> Iterator { return {this, steps()}; }

  ~UniformGrid() = default;
};
