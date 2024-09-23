using Documenter, WaveResolvingBQ

makedocs(
    modules = [WaveResolvingBQ],
    format = :html,
    sitename = "WaveResolvingBQ Documentation",
    pages = [
        "Home" => "index.md",
        "Usage" => "usage.md",
        # Add more pages as needed
    ],
    clean = true
)