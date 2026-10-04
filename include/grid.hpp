#pragma once

#include <cmath>
#include <concepts>
#include <format>
#include <ranges>
#include <stdexcept>
#include <utility>

template <std::floating_point F, std::unsigned_integral U = std::uint32_t>
class UniformGrid final {
private:
    F start_;
    F end_;
    F step_;
    U steps_;

    auto lerp(F pos) const noexcept -> F {
        return std::lerp(this->start_, this->end_, pos / steps_);
    }

    template <std::invocable<U> Fn>
    static auto indexed(U count, Fn fun) {
        return std::views::iota(static_cast<U>(0), count) | std::views::transform(fun);
    }

    template <std::integral I>
    static auto checkIndex(I index, U upper_bound) -> void {
        if (index < 0 or index > upper_bound) {
            throw std::invalid_argument(
                std::format("Index {} is out of bounds [0, {}]", index, upper_bound));
        }
    }

public:
    UniformGrid()                                        = delete;
    UniformGrid(const UniformGrid &)                     = default;
    UniformGrid(UniformGrid &&)                          = default;
    auto operator=(const UniformGrid &) -> UniformGrid & = default;
    auto operator=(UniformGrid &&) -> UniformGrid &      = default;

    template <std::integral I>
        requires std::convertible_to<I, U>
    UniformGrid(F start, F end, I steps)
        : start_(start)
        , end_(end)
        , steps_(steps)
        , step_((end - start) / steps) {
        if (!(start < end)) {
            throw std::invalid_argument("Start must be less than end");
        }
        if (!std::in_range<U>(steps)) {
            throw std::invalid_argument(
                std::format(
                    "Can not convert values {} to type {} without losses",
                    steps,
                    typeid(U).name()));
        }
    }

    [[nodiscard]] auto step() const noexcept -> F {
        return this->step_;
    }

    [[nodiscard]] auto steps() const noexcept -> U {
        return this->steps_;
    }

    [[nodiscard]] auto nodes() const noexcept -> U {
        return this->steps_ + 1;
    }

    template <std::integral I>
    auto region(I index) const -> std::pair<F, F> {
        checkIndex(index, this->steps_);
        return {this->lerp(index), this->lerp(index + 1)};
    }

    template <std::integral I>
    auto node(I index) const -> F {
        checkIndex(index, this->steps_);
        return this->lerp(index);
    }

    template <std::integral I>
    auto supNode(I index) const -> F {
        checkIndex(index, this->steps_ - 1);
        auto main_node = this->node(index);
        // NOLINTNEXTLINE (*-magic-numbers)
        return this->lerp(index + F{0.5});
    }

    template <std::integral I>
    auto supRegion(I index) const -> std::pair<F, F> {
        checkIndex(index, this->steps_ - 2);
        // NOLINTNEXTLINE (*-magic-numbers)
        return {this->lerp(index + F{0.5}), this->lerp(index + F{1.5})};
    }

    auto nodeRange() const {
        return indexed(this->steps_ + 1, [grid = *this](U index) -> F {
            return grid.lerp(index);
        });
    }

    auto supNodeRange() const {
        return indexed(this->steps_, [grid = *this](U index) -> F {
            return grid.lerp(index + F{0.5});  // NOLINT(*-magic-numbers)
        });
    }

    auto regionRange() const {
        return indexed(this->steps_, [grid = *this](U index) -> std::pair<F, F> {
            return {grid.lerp(index), grid.lerp(index + 1)};
        });
    }

    auto supRegionRange() const {
        return indexed(this->steps_ - 1, [grid = *this](U index) -> std::pair<F, F> {
            return {
                grid.lerp(index + F{0.5}),  // NOLINT(*-magic-numbers)
                grid.lerp(index + F{1.5}),  // NOLINT(*-magic-numbers)
            };
        });
    }

    ~UniformGrid() = default;
};
