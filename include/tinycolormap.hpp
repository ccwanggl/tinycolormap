/*
 MIT License

 Copyright (c) 2018-2020 Yuki Koyama

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

 -------------------------------------------------------------------------------

 The lookup table for Turbo is derived by Shot511 in his PR,
 https://github.com/yuki-koyama/tinycolormap/pull/27 , from
 https://gist.github.com/mikhailov-work/6a308c20e494d9e0ccc29036b28faa7a , which
 is released by Anton Mikhailov, copyrighted by Google LLC, and licensed under
 the Apache 2.0 license. To the best of our knowledge, the Apache 2.0 license is
 compatible with the MIT license, and thus we release the merged entire code
 under the MIT license. The license notice for Anton's code is posted here:

 Copyright 2019 Google LLC.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 */

#ifndef TINYCOLORMAP_HPP_
#define TINYCOLORMAP_HPP_

#include <cmath>
#include <cstdint>
#include <algorithm>

#if defined(TINYCOLORMAP_WITH_EIGEN)
#include <Eigen/Core>
#endif

#if defined(TINYCOLORMAP_WITH_QT5)
#include <QColor>
#endif

#if defined(TINYCOLORMAP_WITH_QT5) && defined(TINYCOLORMAP_WITH_EIGEN)
#include <QImage>
#include <QString>
#endif

#if defined(TINYCOLORMAP_WITH_GLM)
#include <glm/vec3.hpp>
#endif

namespace tinycolormap
{
    //////////////////////////////////////////////////////////////////////////////////
    // Interface
    //////////////////////////////////////////////////////////////////////////////////

    enum class ColormapType
    {
        Parula, Heat, Jet, Turbo, Hot, Gray, Magma, BlackBody, Inferno, Plasma, Viridis, Cividis, Github, Cubehelix,HSV
    };

    struct Color
    {
        explicit constexpr Color(double gray) noexcept : data{ gray, gray, gray } {}
        constexpr Color(double r, double g, double b) noexcept : data{ r, g, b } {}

        double data[3];

        double& r() noexcept { return data[0]; }
        double& g() noexcept { return data[1]; }
        double& b() noexcept { return data[2]; }
        constexpr double r() const noexcept { return data[0]; }
        constexpr double g() const noexcept { return data[1]; }
        constexpr double b() const noexcept { return data[2]; }

        constexpr uint8_t ri() const noexcept { return static_cast<uint8_t>(data[0] * 255.0); }
        constexpr uint8_t gi() const noexcept { return static_cast<uint8_t>(data[1] * 255.0); }
        constexpr uint8_t bi() const noexcept { return static_cast<uint8_t>(data[2] * 255.0); }

        double& operator[](std::size_t n) noexcept { return data[n]; }
        constexpr double operator[](std::size_t n) const noexcept { return data[n]; }
        double& operator()(std::size_t n) noexcept { return data[n]; }
        constexpr double operator()(std::size_t n) const noexcept { return data[n]; }

        friend constexpr Color operator+(const Color& c0, const Color& c1) noexcept
        {
            return { c0.r() + c1.r(), c0.g() + c1.g(), c0.b() + c1.b() };
        }

        friend constexpr Color operator*(double s, const Color& c) noexcept
        {
            return { s * c.r(), s * c.g(), s * c.b() };
        }
     
#if defined(TINYCOLORMAP_WITH_QT5)
        QColor ConvertToQColor() const { return QColor(data[0] * 255.0, data[1] * 255.0, data[2] * 255.0); }
#endif
#if defined(TINYCOLORMAP_WITH_EIGEN)
        Eigen::Vector3d ConvertToEigen() const { return Eigen::Vector3d(data[0], data[1], data[2]); }
#endif
#if defined(TINYCOLORMAP_WITH_GLM)
        glm::vec3 ConvertToGLM() const { return glm::vec3(data[0], data[1], data[2]); }
#endif
    };

    inline Color GetColor(double x, ColormapType type = ColormapType::Viridis);
    inline Color GetQuantizedColor(double x, unsigned int num_levels, ColormapType type = ColormapType::Viridis);
    inline Color GetParulaColor(double x);
    inline Color GetHeatColor(double x);
    inline Color GetJetColor(double x);
    inline Color GetTurboColor(double x);
    inline Color GetHotColor(double x);
    inline constexpr Color GetGrayColor(double x) noexcept;
    inline Color GetMagmaColor(double x);
    inline Color GetBlackBodyColor(double x);
    inline Color GetInfernoColor(double x);
    inline Color GetPlasmaColor(double x);
    inline Color GetViridisColor(double x);
    inline Color GetCividisColor(double x);
    inline Color GetGithubColor(double x);
    inline Color GetCubehelixColor(double x);
    inline Color GetHSVColor(double x);

#if defined(TINYCOLORMAP_WITH_QT5) && defined(TINYCOLORMAP_WITH_EIGEN)
    inline QImage CreateMatrixVisualization(const Eigen::MatrixXd& matrix);
    inline void ExportMatrixVisualization(const Eigen::MatrixXd& matrix, const std::string& path);
#endif

    //////////////////////////////////////////////////////////////////////////////////
    // Private Implementation - public usage is not intended
    //////////////////////////////////////////////////////////////////////////////////

    namespace internal
    {
        inline constexpr double Clamp01(double x) noexcept
        {
            return (x < 0.0) ? 0.0 : (x > 1.0) ? 1.0 : x;
        }
        
        // A helper function to calculate linear interpolation
        template <std::size_t N>
        Color CalcLerp(double x, const Color (&data)[N])
        {
            const double a  = Clamp01(x) * (N - 1);
            const double i  = std::floor(a);
            const double t  = a - i;
            const Color& c0 = data[static_cast<std::size_t>(i)];
            const Color& c1 = data[static_cast<std::size_t>(std::ceil(a))];

            return (1.0 - t) * c0 + t * c1;
        }

        inline double QuantizeArgument(double x, unsigned int num_levels)
        {
            // Clamp num_classes to range [1, 255].
            num_levels = (std::max)(1u, (std::min)(num_levels, 255u));

            const double interval_length = 255.0 / num_levels;

            // Calculate index of the interval to which the given x belongs to.
            // Substracting eps prevents getting out of bounds index.
            constexpr double eps = 0.0005;
            const unsigned int index = static_cast<unsigned int>((x * 255.0 - eps) / interval_length);

            // Calculate upper and lower bounds of the given interval.
            const unsigned int upper_boundary = static_cast<unsigned int>(index * interval_length + interval_length);
            const unsigned int lower_boundary = static_cast<unsigned int>(upper_boundary - interval_length);

            // Get middle "coordinate" of the given interval and move it back to [0.0, 1.0] interval.
            const double xx = static_cast<double>(upper_boundary + lower_boundary) * 0.5 / 255.0;

            return xx;
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Public Implementation
    //////////////////////////////////////////////////////////////////////////////////

    inline Color GetColor(double x, ColormapType type)
    {
        switch (type)
        {
            case ColormapType::Parula:
                return GetParulaColor(x);
            case ColormapType::Heat:
                return GetHeatColor(x);
            case ColormapType::Jet:
                return GetJetColor(x);
            case ColormapType::Turbo:
                return GetTurboColor(x);
            case ColormapType::Hot:
                return GetHotColor(x);
            case ColormapType::Gray:
                return GetGrayColor(x);
            case ColormapType::Magma:
                return GetMagmaColor(x);
            case ColormapType::BlackBody:
                return GetBlackBodyColor(x);
            case ColormapType::Inferno:
                return GetInfernoColor(x);
            case ColormapType::Plasma:
                return GetPlasmaColor(x);
            case ColormapType::Viridis:
                return GetViridisColor(x);
            case ColormapType::Cividis:
                return GetCividisColor(x);
            case ColormapType::Github:
                return GetGithubColor(x);
            case ColormapType::Cubehelix:
                return GetCubehelixColor(x);
            case ColormapType::HSV:
                return GetHSVColor(x);
            default:
                break;
        }

        return GetViridisColor(x);
    }

    inline Color GetQuantizedColor(double x, unsigned int num_levels, ColormapType type)
    {
        return GetColor(internal::QuantizeArgument(x, num_levels), type);
    }

    inline Color GetParulaColor(double x)
    {
        constexpr Color data[] =
        {
          #include "ParulaColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetHeatColor(double x)
    {
        constexpr Color data[] =
        {
          #include "HeatColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetJetColor(double x)
    {
        constexpr Color data[] =
        {
          #include "JetColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetTurboColor(double x)
    {
        constexpr Color data[] =
        {
          #include "TurboColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetHotColor(double x)
    {
        x = internal::Clamp01(x);

        constexpr Color r{ 1.0, 0.0, 0.0 };
        constexpr Color g{ 0.0, 1.0, 0.0 };
        constexpr Color b{ 0.0, 0.0, 1.0 };

        if (x < 0.4)
        {
            const double t = x / 0.4;
            return t * r;
        }
        else if (x < 0.8)
        {
            const double t = (x - 0.4) / (0.8 - 0.4);
            return r + t * g;
        }
        else
        {
            const double t = (x - 0.8) / (1.0 - 0.8);
            return r + g + t * b;
        }
    }

    inline constexpr Color GetGrayColor(double x) noexcept
    {
        return Color{ 1.0 - internal::Clamp01(x) };
    }
    inline Color GetBlackBodyColor(double x)
    {
        constexpr Color data[] =
        {
          #include "BlackBodyColor.inc"
		};

			return internal::CalcLerp(x, data);
		}


    inline Color GetMagmaColor(double x)
    {
        constexpr Color data[] =
        {
          #include "MagmaColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetInfernoColor(double x)
    {
        constexpr Color data[] =
        {
          #include "InfernoColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetPlasmaColor(double x)
    {
        constexpr Color data[] =
        {
          #include "PlasmaColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetViridisColor(double x)
    {
        constexpr Color data[] =
        {
          #include "ViridisColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetCividisColor(double x)
    {
        constexpr Color data[] =
        {
          #include "CividisColor.inc"
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetGithubColor(double x)
    {
        constexpr Color data[] =
        {
            { 0.933333, 0.933333, 0.933333 },
            { 0.776470, 0.894117, 0.545098 },
            { 0.482352, 0.788235, 0.435294 },
            { 0.137254, 0.603921, 0.231372 },
            { 0.098039, 0.380392, 0.152941 }
        };

        return internal::CalcLerp(x, data);
    }

    inline Color GetCubehelixColor(double x)
    {
        constexpr Color data[] =
        {
          #include "CubehelixColor.inc"
        };

        return internal::CalcLerp(x, data);
    }
    inline Color GetHSVColor(double x)
    {
        constexpr Color data[] =
        {
          #include "HSVColor.inc"
		};
			return internal::CalcLerp(x, data);
		}

#if defined(TINYCOLORMAP_WITH_QT5) && defined(TINYCOLORMAP_WITH_EIGEN)
    inline QImage CreateMatrixVisualization(const Eigen::MatrixXd& matrix)
    {
        const int w = matrix.cols();
        const int h = matrix.rows();
        const double max_coeff = matrix.maxCoeff();
        const double min_coeff = matrix.minCoeff();
        const Eigen::MatrixXd normalized = (1.0 / (max_coeff - min_coeff)) * (matrix - Eigen::MatrixXd::Constant(h, w, min_coeff));

        QImage image(w, h, QImage::Format_ARGB32);
        for (int x = 0; x < w; ++ x)
        {
            for (int y = 0; y < h; ++ y)
            {
                const QColor color = tinycolormap::GetColor(normalized(y, x)).ConvertToQColor();
                image.setPixel(x, y, color.rgb());
            }
        }

        return image;
    }

    inline void ExportMatrixVisualization(const Eigen::MatrixXd& matrix, const std::string& path)
    {
        CreateMatrixVisualization(matrix).save(QString::fromStdString(path));
    }
#endif
}

#endif
