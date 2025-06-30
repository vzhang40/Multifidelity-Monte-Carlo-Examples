function plotModels(models, a, b, fig)
    figure(fig); clf(fig);
    x = linspace(a, b, 100);
    legendtext = cell(size({models.f}));
    legendtext{1} = "Truth $f^{(1)}$";
    plot(x, models(1).f(x))
    hold on
    for i = 2:length({models.f})
        plot(x, models(i).f(x))
        legendtext{i} = "$f^{(" + i + ")}$";
    end
    legend(legendtext, "Interpreter", "latex")
    xlabel("$x$", "Interpreter", "latex")
    ylabel("$f^{(i)}(x)$", "Interpreter", "latex")
    title("Truth and Surrogate Models")
end
